"""The aligned-trait task: both input formats, allele flips, and the upstream-bug guards."""

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.fast_ct_ldsc.aligned_trait_task import (
    AlignedTraitArrays,
    FastCtLdscAlignedTraitTask,
    read_aligned_trait,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.alignment_guards import AlignmentGuards
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
    SumstatsSource,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
)
from test_mecfs_bio.unit.build_system.task.fast_ct_ldsc.fast_ct_ldsc_fixtures import (
    SyntheticIndex,
    build_synthetic_index,
    parquet_trait_meta,
    pickle_trait_meta,
    register_frame,
    register_pickle,
    run_task,
    trait_frame,
    trait_sumstats,
)

_N_SNPS = 400
_N = 10_000.0


def _truth_z(n_snps: int = _N_SNPS) -> np.ndarray:
    return np.random.default_rng(1).standard_normal(n_snps)


def _run_parquet(
    tmp_path: Path,
    frame: pl.DataFrame,
    *,
    genome_build: GenomeBuild | None = "19",
    guards: AlignmentGuards = AlignmentGuards(),
    source_guards: AlignmentGuards | None = None,
) -> AlignedTraitArrays:
    idx = build_synthetic_index(tmp_path, _N_SNPS, seed=0)
    trait_task = FakeTask(
        meta=parquet_trait_meta(asset_id="trait", genome_build=genome_build)
    )
    register_frame(tmp_path, idx.assets, trait_task, frame)
    return _align(tmp_path, idx, trait_task, guards=guards, source_guards=source_guards)


def _run_pickle(
    tmp_path: Path,
    frame: pl.DataFrame,
    *,
    meta_build: GenomeBuild | None,
    pickle_build: GenomeBuild | None,
) -> AlignedTraitArrays:
    idx = build_synthetic_index(tmp_path, _N_SNPS, seed=0)
    trait_task = FakeTask(
        meta=pickle_trait_meta(asset_id="trait", genome_build=meta_build)
    )
    register_pickle(
        tmp_path, idx.assets, trait_task, trait_sumstats(frame, pickle_build)
    )
    return _align(
        tmp_path, idx, trait_task, guards=AlignmentGuards(), source_guards=None
    )


def _align(
    tmp_path: Path,
    idx: SyntheticIndex,
    trait_task: FakeTask,
    *,
    guards: AlignmentGuards,
    source_guards: AlignmentGuards | None,
) -> AlignedTraitArrays:
    source = SumstatsSource(
        task=trait_task,
        alias="T",
        sample_info=QuantPhenotype(),
        alignment_guards=source_guards,
    )
    task = FastCtLdscAlignedTraitTask.create(
        "aligned_t", source=source, index_task=idx.index_task, guards=guards
    )
    idx.assets[task.asset_id] = run_task(tmp_path, task, idx.assets)
    return read_aligned_trait(task, lambda asset_id: idx.assets[asset_id])


def _variants(tmp_path: Path) -> pl.DataFrame:
    return build_synthetic_index(tmp_path / "v", _N_SNPS, seed=0).variants


@pytest.fixture
def variants(tmp_path: Path) -> pl.DataFrame:
    (tmp_path / "v").mkdir()
    return _variants(tmp_path)


def test_parquet_trait_is_aligned_with_flips_and_writes_diagnostics(
    tmp_path: Path, variants: pl.DataFrame
) -> None:
    z = _truth_z()
    arrays = _run_parquet(tmp_path, trait_frame(variants, z, n=_N))
    np.testing.assert_allclose(arrays.z, z)
    np.testing.assert_allclose(arrays.n, _N)
    assert arrays.diagnostics.n_matched == _N_SNPS
    assert arrays.diagnostics.position_concordance == 1.0


def test_gwaslab_pickle_matches_parquet_input(
    tmp_path: Path, variants: pl.DataFrame
) -> None:
    frame = trait_frame(variants, _truth_z(), n=_N)
    (tmp_path / "p").mkdir()
    (tmp_path / "q").mkdir()
    from_parquet = _run_parquet(tmp_path / "p", frame)
    from_pickle = _run_pickle(tmp_path / "q", frame, meta_build="19", pickle_build="19")
    np.testing.assert_array_equal(from_pickle.z, from_parquet.z)
    np.testing.assert_array_equal(from_pickle.n, from_parquet.n)


def test_shuffled_rsids_trip_mismatch_guard(
    tmp_path: Path, variants: pl.DataFrame
) -> None:
    frame = trait_frame(variants, _truth_z(), n=_N)
    shuffled = frame.with_columns(
        pl.Series(
            GWASLAB_RSID_COL,
            np.random.default_rng(2).permutation(frame[GWASLAB_RSID_COL].to_numpy()),
        )
    )
    with pytest.raises(AssertionError):
        _run_parquet(
            tmp_path,
            shuffled,
            genome_build=None,
            guards=AlignmentGuards(min_coverage=0.0),
        )


def test_sparse_trait_trips_coverage_guard(
    tmp_path: Path, variants: pl.DataFrame
) -> None:
    z = _truth_z()
    z[np.arange(_N_SNPS) % 10 != 0] = np.nan
    with pytest.raises(AssertionError):
        _run_parquet(tmp_path, trait_frame(variants, z, n=_N))


def _shifted(variants: pl.DataFrame) -> pl.DataFrame:
    frame = trait_frame(variants, _truth_z(), n=_N)
    return frame.with_columns(pl.col(GWASLAB_POS_COL) + 7)


def test_shifted_positions_trip_concordance_guard_when_declared_build_19(
    tmp_path: Path, variants: pl.DataFrame
) -> None:
    with pytest.raises(AssertionError):
        _run_parquet(tmp_path, _shifted(variants), genome_build="19")


def test_gwaslab_pickle_without_meta_build_uses_pickle_build(
    tmp_path: Path, variants: pl.DataFrame
) -> None:
    with pytest.raises(AssertionError):
        _run_pickle(tmp_path, _shifted(variants), meta_build=None, pickle_build="19")


def test_gwaslab_meta_build_contradicting_pickle_build_fails(
    tmp_path: Path, variants: pl.DataFrame
) -> None:
    frame = trait_frame(variants, _truth_z(), n=_N)
    with pytest.raises(AssertionError):
        _run_pickle(tmp_path, frame, meta_build="38", pickle_build="19")


@pytest.mark.parametrize("genome_build", [None, "38"])
def test_shifted_positions_pass_without_declared_build_19(
    tmp_path: Path, variants: pl.DataFrame, genome_build: GenomeBuild | None
) -> None:
    arrays = _run_parquet(tmp_path, _shifted(variants), genome_build=genome_build)
    assert arrays.diagnostics.position_concordance == 0.0


def test_per_source_override_lets_a_tripping_trait_through(
    tmp_path: Path, variants: pl.DataFrame
) -> None:
    z = _truth_z()
    z[np.arange(_N_SNPS) % 10 != 0] = np.nan
    arrays = _run_parquet(
        tmp_path,
        trait_frame(variants, z, n=_N),
        source_guards=AlignmentGuards(min_coverage=0.05),
    )
    assert arrays.diagnostics.n_matched == _N_SNPS // 10


def test_null_n_with_finite_z_fails(tmp_path: Path, variants: pl.DataFrame) -> None:
    frame = trait_frame(variants, _truth_z(), n=_N).with_columns(
        pl.when(pl.int_range(pl.len()) == 3)
        .then(None)
        .otherwise(pl.col(GWASLAB_SAMPLE_SIZE_COLUMN))
        .alias(GWASLAB_SAMPLE_SIZE_COLUMN)
    )
    with pytest.raises(AssertionError):
        _run_parquet(tmp_path, frame)


def test_aligned_rsids_written_in_context_order(
    tmp_path: Path, variants: pl.DataFrame
) -> None:
    # Reverse the input rows: the output must still follow the context order.
    frame = trait_frame(variants, _truth_z(), n=_N).reverse()
    arrays = _run_parquet(tmp_path, frame)
    assert arrays.rsid.to_list() == variants[GWASLAB_RSID_COL].to_list()
