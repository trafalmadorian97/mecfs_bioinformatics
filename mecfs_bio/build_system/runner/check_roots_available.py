from typing import Sequence

from mecfs_bio.build_system.rebuilder.metadata_to_path.remapping_meta_to_path import (
    PathRemapRule,
)


class RemapRootUnavailableError(RuntimeError):
    """
    A configured remap root is not present, so the subtrees routed there can be neither
    read nor written.
    """


def _unavailable_root_message(rule: PathRemapRule) -> str:
    if rule.root.exists():
        cause = f"{rule.root} exists but is not a directory."
    else:
        cause = f"{rule.root} does not exist. "
    prefixes = ", ".join(str(prefix) for prefix in rule.prefixes)
    return (
        f"Asset store remap root unavailable: {cause}\n"
        f"Subtrees routed there: {prefixes}.\n"
    )


def check_remap_roots_available(rules: Sequence[PathRemapRule]) -> None:
    """
    Fail fast if any configured remap root is missing.

    """
    unavailable = [rule for rule in rules if not rule.root.is_dir()]
    if not unavailable:
        return
    raise RemapRootUnavailableError(
        "\n\n".join(_unavailable_root_message(rule) for rule in unavailable)
    )
