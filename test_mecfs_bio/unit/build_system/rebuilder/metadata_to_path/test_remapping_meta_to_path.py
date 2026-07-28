from pathlib import Path, PurePath

import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.metadata_to_path.remapping_meta_to_path import (
    MIGRATE_COMMAND,
    PathRemapRule,
    RemappingMetaToPath,
    RemapRootUnavailableError,
    check_remap_roots_available,
    stale_default_root_dirs,
)
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    SimpleMetaToPath,
    simple_meta_to_relative_path,
)

DEFAULT_ROOT = Path("/default_root")
REMAP_ROOT = Path("/remap_root")

DB_SNP_PREFIX = PurePath("reference_data/db_snp_reference_data")


def _reference_file_meta(group: str, asset_id: str = "some_file") -> ReferenceFileMeta:
    return ReferenceFileMeta(
        group=group,
        sub_group="some_sub_group",
        sub_folder=PurePath("some_sub_folder"),
        extension=".vcf.gz",
        id=AssetId(asset_id),
    )


ALL_META_KINDS = [
    SimpleFileMeta(AssetId("a_file")),
    SimpleDirectoryMeta(AssetId("a_directory")),
    _reference_file_meta("db_snp_reference_data"),
    _reference_file_meta("linkage_disequilibrium_scores"),
]


@pytest.mark.parametrize(argnames="meta", argvalues=ALL_META_KINDS)
def test_no_rules_is_equivalent_to_simple_meta_to_path(meta: Meta) -> None:
    """
    Remapping with an empty rule set must not perturb the existing store layout.
    """
    remapping = RemappingMetaToPath(default_root=DEFAULT_ROOT, rules=())
    assert remapping(meta) == SimpleMetaToPath(root=DEFAULT_ROOT)(meta)


def test_matched_prefix_is_routed_to_the_rule_root() -> None:
    meta = _reference_file_meta("db_snp_reference_data")
    remapping = RemappingMetaToPath(
        default_root=DEFAULT_ROOT,
        rules=(PathRemapRule(root=REMAP_ROOT, prefixes=(DB_SNP_PREFIX,)),),
    )
    assert remapping(meta) == REMAP_ROOT / simple_meta_to_relative_path(meta)


def test_unmatched_prefix_stays_under_the_default_root() -> None:
    meta = _reference_file_meta("linkage_disequilibrium_scores")
    remapping = RemappingMetaToPath(
        default_root=DEFAULT_ROOT,
        rules=(PathRemapRule(root=REMAP_ROOT, prefixes=(DB_SNP_PREFIX,)),),
    )
    assert remapping(meta) == SimpleMetaToPath(root=DEFAULT_ROOT)(meta)


def test_prefix_matching_is_component_wise() -> None:
    """
    A string prefix match would wrongly capture a sibling whose name merely starts with
    the prefix, sending it to a different filesystem than intended.
    """
    sibling = _reference_file_meta("db_snp_reference_data_extra")
    remapping = RemappingMetaToPath(
        default_root=DEFAULT_ROOT,
        rules=(PathRemapRule(root=REMAP_ROOT, prefixes=(DB_SNP_PREFIX,)),),
    )
    assert remapping(sibling) == SimpleMetaToPath(root=DEFAULT_ROOT)(sibling)


def test_overlapping_prefixes_across_rules_are_rejected() -> None:
    with pytest.raises(AssertionError):
        RemappingMetaToPath(
            default_root=DEFAULT_ROOT,
            rules=(
                PathRemapRule(root=REMAP_ROOT, prefixes=(PurePath("reference_data"),)),
                PathRemapRule(root=Path("/other_root"), prefixes=(DB_SNP_PREFIX,)),
            ),
        )


def test_overlapping_prefixes_within_a_rule_are_rejected() -> None:
    with pytest.raises(AssertionError):
        RemappingMetaToPath(
            default_root=DEFAULT_ROOT,
            rules=(
                PathRemapRule(
                    root=REMAP_ROOT,
                    prefixes=(PurePath("reference_data"), DB_SNP_PREFIX),
                ),
            ),
        )


def test_absolute_prefix_is_rejected() -> None:
    with pytest.raises(AssertionError):
        PathRemapRule(root=REMAP_ROOT, prefixes=(PurePath("/reference_data"),))


def test_empty_prefix_list_is_rejected() -> None:
    with pytest.raises(AssertionError):
        PathRemapRule(root=REMAP_ROOT, prefixes=())


def test_tuple_from_config_parses_a_realistic_section() -> None:
    rules = PathRemapRule.tuple_from_config(
        [
            {
                "root": "/mnt/d/asset_store_remote",
                "prefixes": [
                    "reference_data/db_snp_reference_data",
                    "reference_data/genome_annotations",
                ],
            }
        ]
    )
    assert rules == (
        PathRemapRule(
            root=Path("/mnt/d/asset_store_remote"),
            prefixes=(
                PurePath("reference_data/db_snp_reference_data"),
                PurePath("reference_data/genome_annotations"),
            ),
        ),
    )


def test_tuple_from_config_rejects_unknown_keys() -> None:
    """
    A typo in a machine-local config should fail loudly rather than silently disabling a
    rule and sending assets back to the disk the remapping was meant to relieve.
    """
    with pytest.raises(AssertionError):
        PathRemapRule.tuple_from_config(
            [{"root": "/mnt/d", "prefix": ["reference_data/db_snp_reference_data"]}]
        )


def test_available_remap_root_passes_the_check(tmp_path: Path) -> None:
    check_remap_roots_available(
        (PathRemapRule(root=tmp_path, prefixes=(DB_SNP_PREFIX,)),)
    )


def test_no_rules_needs_no_roots() -> None:
    check_remap_roots_available(())


def test_missing_remap_root_is_reported_with_enough_detail_to_act_on(
    tmp_path: Path,
) -> None:
    """
    This is the detached-drive case, which is otherwise silent: the assets merely look
    unbuilt.  The message has to name the root, say what is routed there, and give the
    reader somewhere to go.
    """
    missing_root = tmp_path / "not_mounted"
    with pytest.raises(RemapRootUnavailableError) as raised:
        check_remap_roots_available(
            (PathRemapRule(root=missing_root, prefixes=(DB_SNP_PREFIX,)),)
        )

    message = str(raised.value)
    assert str(missing_root) in message
    assert str(DB_SNP_PREFIX) in message
    assert MIGRATE_COMMAND in message


def test_remap_root_that_is_a_file_is_rejected(tmp_path: Path) -> None:
    not_a_directory = tmp_path / "a_file"
    not_a_directory.write_text("")
    with pytest.raises(RemapRootUnavailableError):
        check_remap_roots_available(
            (PathRemapRule(root=not_a_directory, prefixes=(DB_SNP_PREFIX,)),)
        )


def test_every_missing_root_is_reported_at_once(tmp_path: Path) -> None:
    """
    Reporting one root at a time would have the user attach a drive, restart, and hit the
    next failure.
    """
    first = tmp_path / "first_missing"
    second = tmp_path / "second_missing"
    with pytest.raises(RemapRootUnavailableError) as raised:
        check_remap_roots_available(
            (
                PathRemapRule(root=first, prefixes=(DB_SNP_PREFIX,)),
                PathRemapRule(root=second, prefixes=(PurePath("gwas/some_trait"),)),
            )
        )

    message = str(raised.value)
    assert str(first) in message
    assert str(second) in message


def test_stale_default_root_dirs_reports_unmigrated_subtrees(tmp_path: Path) -> None:
    rules = (
        PathRemapRule(
            root=REMAP_ROOT,
            prefixes=(DB_SNP_PREFIX, PurePath("reference_data/genome_annotations")),
        ),
    )
    assert stale_default_root_dirs(tmp_path, rules) == []
    (tmp_path / DB_SNP_PREFIX).mkdir(parents=True)
    assert stale_default_root_dirs(tmp_path, rules) == [tmp_path / DB_SNP_PREFIX]
