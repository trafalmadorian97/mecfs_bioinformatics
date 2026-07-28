"""
Spread one logical asset store across several filesystems.

Motivation: the asset store outgrows the disk it lives on long before the machine runs
out of storage in total.  Relocating the whole store to a large external drive is a bad
trade, because such drives are usually reached over a mount whose dominant cost is
per-file latency rather than bulk throughput.  Subtrees made of thousands of small files
become painfully slow there, while a subtree of a few enormous files is barely affected.

The remedy is to route individual subtrees, not the whole store.  Because
simple_meta_to_relative_path is root-free, an asset's location is a root plus a fixed
relative path, so routing is entirely a question of which root to prepend.

Choosing what to remap: prefer subtrees that are large, made of FEW files, and read
infrequently.  Measured examples from this repo's own store show how far apart the two
extremes sit:

    reference_data/db_snp_reference_data      85G in 9 files       ideal to remap
    reference_data/genome_annotations         21G in 2 files       ideal to remap
    reference_data/linkage_disequilibrium_scores
                                              51G in 186,105 files must stay local

Two consequences worth knowing before remapping a subtree:

    Finalizing a remapped asset is no longer an instant os.rename.  sandboxed_execute
    builds into a scratch directory on the local disk, so moving the result to a remapped
    root crosses a filesystem boundary and falls back to _cross_device_move, which copies
    every byte.  This is correct but slow, and is a further reason to remap only subtrees
    that are rarely rebuilt.

    Remapping lowers steady-state usage of the default root, not peak usage.  The scratch
    directory still lives on the local disk, so building a large remapped asset still
    needs transient room for it there.  TMPDIR is the separate knob for that.

Relocating an asset does not invalidate the build cache: verifying traces are content
hashes and carry no location.  Moving the config without moving the files, however, makes
existing assets invisible to get_asset_if_exists and forces a rebuild, so a rule change
must be paired with an actual migration of the data.
"""

from pathlib import Path, PurePath
from typing import Any, Callable, Mapping, Sequence

import structlog
from attrs import frozen

from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.metadata_to_path.base_meta_to_path import (
    MetaToPath,
)
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    simple_meta_to_relative_path,
)

logger = structlog.get_logger(__name__)

_ROOT_KEY = "root"
_PREFIXES_KEY = "prefixes"
_ALLOWED_RULE_KEYS = frozenset({_ROOT_KEY, _PREFIXES_KEY})

# Repo-level command that reconciles the store on disk with the configured rules.
MIGRATE_COMMAND = "pixi r invoke migrate-asset-store"


@frozen
class PathRemapRule:
    """
    Send every asset whose store-relative path starts with one of prefixes to root
    instead of to the default store root.

    In yaml, as consumed from default_runner_config.yaml:

        path_remap:
          - root: /mnt/d/asset_store_remote
            prefixes:
              - reference_data/db_snp_reference_data
              - reference_data/genome_annotations

    A prefix is matched component-wise, never as a string prefix, so a rule on
    reference_data/db_snp does not also capture reference_data/db_snp_extra.

    Prefixes are ordinary store-relative paths, so any level of the layout can be
    selected.  gwas/<trait> works just as well as reference_data/<group>, which makes it
    possible to archive a trait that is no longer under active analysis.
    """

    root: Path
    prefixes: tuple[PurePath, ...]

    def __attrs_post_init__(self) -> None:
        assert self.prefixes, f"Remap rule for root {self.root} lists no prefixes."
        for prefix in self.prefixes:
            assert prefix.parts, (
                f"Empty remap prefix for root {self.root}; a prefix must select a subtree."
            )
            assert not prefix.is_absolute(), (
                f"Remap prefix {prefix} must be relative to the asset store root, "
                f"not absolute."
            )

    def matches(self, relative: PurePath) -> bool:
        return any(
            relative.parts[: len(prefix.parts)] == prefix.parts
            for prefix in self.prefixes
        )

    @classmethod
    def tuple_from_config(
        cls, entries: Sequence[Mapping[str, Any]]
    ) -> tuple["PathRemapRule", ...]:
        """
        Build rules from the path_remap section of a parsed yaml config.  Unrecognized
        keys are an error rather than being ignored, so a typo in a machine-local config
        surfaces immediately instead of silently disabling a rule.
        """
        rules = []
        for entry in entries:
            unknown = set(entry) - _ALLOWED_RULE_KEYS
            assert not unknown, (
                f"Unknown key(s) {sorted(unknown)} in a path_remap entry; "
                f"expected only {sorted(_ALLOWED_RULE_KEYS)}."
            )
            assert _ROOT_KEY in entry, f"path_remap entry {entry} is missing a root."
            assert _PREFIXES_KEY in entry, (
                f"path_remap entry {entry} is missing prefixes."
            )
            rules.append(
                cls(
                    root=Path(entry[_ROOT_KEY]),
                    prefixes=tuple(PurePath(p) for p in entry[_PREFIXES_KEY]),
                )
            )
        return tuple(rules)


@frozen
class RemappingMetaToPath(MetaToPath):
    """
    Place assets under default_root, except for those selected by a rule in rules, which
    go under that rule's root instead.

    With no rules this is exactly equivalent to SimpleMetaToPath(root=default_root).
    """

    default_root: Path
    rules: tuple[PathRemapRule, ...]
    relative_path: Callable[[Meta], PurePath] = simple_meta_to_relative_path

    def __attrs_post_init__(self) -> None:
        _check_prefixes_disjoint(self.rules)

    def __call__(self, m: Meta) -> Path:
        relative = self.relative_path(m)
        for rule in self.rules:
            if rule.matches(relative):
                return rule.root / relative
        return self.default_root / relative


def _check_prefixes_disjoint(rules: Sequence[PathRemapRule]) -> None:
    """
    Reject rule sets in which one prefix is a prefix of another, whether within a rule or
    across rules.  Overlapping prefixes would make an asset's location depend on rule
    order, which is not something a config file should have to express.
    """
    seen: list[tuple[PurePath, Path]] = [
        (prefix, rule.root) for rule in rules for prefix in rule.prefixes
    ]
    for i, (prefix, root) in enumerate(seen):
        for other_prefix, other_root in seen[i + 1 :]:
            shorter, longer = sorted((prefix, other_prefix), key=lambda p: len(p.parts))
            assert longer.parts[: len(shorter.parts)] != shorter.parts, (
                f"Remap prefixes overlap: {prefix} (root {root}) and "
                f"{other_prefix} (root {other_root}).  One asset would match both, so "
                f"its location would depend on rule order."
            )


def stale_default_root_dirs(
    default_root: Path, rules: Sequence[PathRemapRule]
) -> list[Path]:
    """
    Return the directories that are still present under default_root even though a rule
    now routes them elsewhere.  A non-empty result means assets have not been migrated
    and will be rebuilt rather than found.
    """
    return [
        default_root / prefix
        for rule in rules
        for prefix in rule.prefixes
        if (default_root / prefix).exists()
    ]


def warn_on_unmigrated_assets(
    default_root: Path, rules: Sequence[PathRemapRule]
) -> None:
    """
    Warn if a subtree that is now routed elsewhere still holds data under default_root.

    That state means the config was changed without moving the data, so those assets are
    invisible to the build system and will be rebuilt from scratch at the new location
    while the old copy keeps occupying the disk the remapping was meant to relieve.

    This is deliberately only a warning: aborting would strand a pipeline that may already
    be hours into a run, and rebuilding is wasteful rather than wrong.  It costs one stat
    per configured prefix, and nothing at all when no rules are configured.
    """
    for stale in stale_default_root_dirs(default_root, rules):
        logger.warning(
            f"Asset store subtree {stale} is still present under the default asset root "
            f"even though path_remap now routes it elsewhere.  Those assets will be "
            f"rebuilt rather than found.  Run '{MIGRATE_COMMAND}' to move them."
        )
