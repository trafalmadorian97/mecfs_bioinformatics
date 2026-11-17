from pathlib import Path, PurePath

import pytest

from mecfs_bio.build_system.rebuilder.metadata_to_path.remap_globs import (
    GlobExpansion,
    RemapGlob,
    StoreMove,
    expand_glob,
    globs_from_config,
    plan_glob_reconcile,
    read_glob_expansions,
    rules_with_expansions,
    write_glob_expansions,
)
from mecfs_bio.build_system.rebuilder.metadata_to_path.remapping_meta_to_path import (
    PathRemapRule,
    RemappingMetaToPath,
)

RAW_GLOB = "gwas/**/raw"
REMAP_ROOT = Path("/remap_root")
OTHER_ROOT = Path("/other_root")
ALS_RAW = PurePath("gwas/als/van_rheenen/raw")
LDL_RAW = PurePath("gwas/ldl/willer/raw")


def _make_dirs(root: Path, *relative: PurePath) -> None:
    for r in relative:
        (root / r).mkdir(parents=True)


def test_expand_glob_returns_outermost_matching_directories(tmp_path: Path) -> None:
    _make_dirs(
        tmp_path,
        ALS_RAW,
        LDL_RAW / "nested" / "raw",
        PurePath("gwas/ldl/willer/processed"),
        PurePath("reference_data/x/raw"),
    )
    (tmp_path / "gwas/ldl/raw").write_text("a file, not a directory")
    assert expand_glob(tmp_path, RemapGlob(root=REMAP_ROOT, pattern=RAW_GLOB)) == (
        ALS_RAW,
        LDL_RAW,
    )


def test_expand_glob_with_absent_subtree_is_empty(tmp_path: Path) -> None:
    assert expand_glob(tmp_path, RemapGlob(root=REMAP_ROOT, pattern=RAW_GLOB)) == ()


@pytest.mark.parametrize(argnames="pattern", argvalues=["**/raw", "/gwas/**/raw", ""])
def test_glob_must_be_relative_with_a_literal_head(pattern: str) -> None:
    with pytest.raises(AssertionError):
        RemapGlob(root=REMAP_ROOT, pattern=pattern)


def test_globs_only_entry_parses_to_globs_and_no_literal_rule() -> None:
    entries = [{"root": str(REMAP_ROOT), "globs": [RAW_GLOB]}]
    assert PathRemapRule.tuple_from_config(entries) == ()
    assert globs_from_config(entries) == (RemapGlob(root=REMAP_ROOT, pattern=RAW_GLOB),)


def test_entry_with_neither_prefixes_nor_globs_is_rejected() -> None:
    with pytest.raises(AssertionError):
        PathRemapRule.tuple_from_config([{"root": str(REMAP_ROOT)}])


def test_expansions_file_round_trips(tmp_path: Path) -> None:
    expansions = (
        GlobExpansion(root=REMAP_ROOT, pattern=RAW_GLOB, prefixes=(ALS_RAW, LDL_RAW)),
    )
    path = tmp_path / "expansions.yaml"
    write_glob_expansions(path, expansions)
    assert read_glob_expansions(path) == expansions


def test_missing_expansions_file_reads_as_empty(tmp_path: Path) -> None:
    assert read_glob_expansions(tmp_path / "absent.yaml") == ()


def test_rules_with_expansions_routes_recorded_prefixes() -> None:
    """
    Includes an expansion whose root has no configured rule: its data still lives there
    until reconcile moves it, so the runner must keep finding it.
    """
    literal = PurePath("reference_data/db_snp")
    rules = rules_with_expansions(
        literal_rules=(PathRemapRule(root=REMAP_ROOT, prefixes=(literal,)),),
        expansions=(
            GlobExpansion(root=REMAP_ROOT, pattern=RAW_GLOB, prefixes=(ALS_RAW,)),
            GlobExpansion(root=OTHER_ROOT, pattern="gwas/ldl/*", prefixes=(LDL_RAW,)),
        ),
    )
    mapper = RemappingMetaToPath(default_root=Path("/default"), rules=rules)
    resolved = {rule.root: set(rule.prefixes) for rule in mapper.rules}
    assert resolved == {REMAP_ROOT: {literal, ALS_RAW}, OTHER_ROOT: {LDL_RAW}}


def test_new_local_match_is_moved_to_the_glob_root(tmp_path: Path) -> None:
    _make_dirs(tmp_path, ALS_RAW)
    plan = plan_glob_reconcile(
        asset_root=tmp_path,
        literal_rules=(),
        globs=(RemapGlob(root=REMAP_ROOT, pattern=RAW_GLOB),),
        previous=(),
    )
    assert plan.moves == (
        StoreMove(prefix=ALS_RAW, source_root=tmp_path, dest_root=REMAP_ROOT),
    )
    assert plan.expansions == (
        GlobExpansion(root=REMAP_ROOT, pattern=RAW_GLOB, prefixes=(ALS_RAW,)),
    )


def test_already_remapped_prefix_is_kept_without_a_move(tmp_path: Path) -> None:
    """
    Remapped data is absent from the default root, so expansion cannot find it there.
    It must still be recorded, or the next reconcile would move it back.
    """
    previous = (GlobExpansion(root=REMAP_ROOT, pattern=RAW_GLOB, prefixes=(ALS_RAW,)),)
    plan = plan_glob_reconcile(
        asset_root=tmp_path,
        literal_rules=(),
        globs=(RemapGlob(root=REMAP_ROOT, pattern=RAW_GLOB),),
        previous=previous,
    )
    assert plan.moves == ()
    assert plan.expansions == previous


def test_removed_glob_moves_data_back_to_the_default_root(tmp_path: Path) -> None:
    plan = plan_glob_reconcile(
        asset_root=tmp_path,
        literal_rules=(),
        globs=(),
        previous=(
            GlobExpansion(root=REMAP_ROOT, pattern=RAW_GLOB, prefixes=(ALS_RAW,)),
        ),
    )
    assert plan.moves == (
        StoreMove(prefix=ALS_RAW, source_root=REMAP_ROOT, dest_root=tmp_path),
    )
    assert plan.expansions == ()


def test_narrowed_glob_moves_back_only_the_dropped_prefixes(tmp_path: Path) -> None:
    plan = plan_glob_reconcile(
        asset_root=tmp_path,
        literal_rules=(),
        globs=(RemapGlob(root=REMAP_ROOT, pattern="gwas/als/**/raw"),),
        previous=(
            GlobExpansion(
                root=REMAP_ROOT, pattern=RAW_GLOB, prefixes=(ALS_RAW, LDL_RAW)
            ),
        ),
    )
    assert plan.moves == (
        StoreMove(prefix=LDL_RAW, source_root=REMAP_ROOT, dest_root=tmp_path),
    )


def test_glob_moved_to_another_root_moves_directly(tmp_path: Path) -> None:
    plan = plan_glob_reconcile(
        asset_root=tmp_path,
        literal_rules=(),
        globs=(RemapGlob(root=OTHER_ROOT, pattern=RAW_GLOB),),
        previous=(
            GlobExpansion(root=REMAP_ROOT, pattern=RAW_GLOB, prefixes=(ALS_RAW,)),
        ),
    )
    assert plan.moves == (
        StoreMove(prefix=ALS_RAW, source_root=REMAP_ROOT, dest_root=OTHER_ROOT),
    )


def test_match_covered_by_same_root_literal_prefix_is_not_recorded(
    tmp_path: Path,
) -> None:
    _make_dirs(tmp_path, ALS_RAW)
    plan = plan_glob_reconcile(
        asset_root=tmp_path,
        literal_rules=(
            PathRemapRule(root=REMAP_ROOT, prefixes=(PurePath("gwas/als"),)),
        ),
        globs=(RemapGlob(root=REMAP_ROOT, pattern=RAW_GLOB),),
        previous=(),
    )
    assert plan.expansions == ()


@pytest.mark.parametrize(
    argnames="literal",
    argvalues=[
        PathRemapRule(root=OTHER_ROOT, prefixes=(PurePath("gwas/als"),)),
        PathRemapRule(root=REMAP_ROOT, prefixes=(ALS_RAW / "inner",)),
    ],
)
def test_match_overlapping_a_literal_prefix_is_rejected(
    tmp_path: Path, literal: PathRemapRule
) -> None:
    _make_dirs(tmp_path, ALS_RAW)
    with pytest.raises(AssertionError):
        plan_glob_reconcile(
            asset_root=tmp_path,
            literal_rules=(literal,),
            globs=(RemapGlob(root=REMAP_ROOT, pattern=RAW_GLOB),),
            previous=(),
        )


def test_globs_selecting_one_subtree_for_different_roots_are_rejected(
    tmp_path: Path,
) -> None:
    _make_dirs(tmp_path, ALS_RAW)
    with pytest.raises(AssertionError):
        plan_glob_reconcile(
            asset_root=tmp_path,
            literal_rules=(),
            globs=(
                RemapGlob(root=REMAP_ROOT, pattern=RAW_GLOB),
                RemapGlob(root=OTHER_ROOT, pattern="gwas/als/*"),
            ),
            previous=(),
        )


def test_failed_moves_stay_recorded_where_their_data_is(tmp_path: Path) -> None:
    """
    A failed reverse move leaves the data on the remap root, so the old record must
    survive.  A failed forward move leaves it on the default root, so no record may be
    written for it.
    """
    _make_dirs(tmp_path, LDL_RAW)
    plan = plan_glob_reconcile(
        asset_root=tmp_path,
        literal_rules=(),
        globs=(RemapGlob(root=REMAP_ROOT, pattern="gwas/ldl/**/raw"),),
        previous=(
            GlobExpansion(root=REMAP_ROOT, pattern=RAW_GLOB, prefixes=(ALS_RAW,)),
        ),
    )
    assert plan.expansions_after(failed=(ALS_RAW, LDL_RAW)) == plan.previous
