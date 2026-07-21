"""Smoke test for PppProteinCrossTraitRgTask over synthetic inputs.

This exercises the task's wiring rather than its arithmetic: the LDSC context, trait alignment,
sample-size lookup, per-protein batching and output assembly all run end to end on dataframes
small enough to write by hand. The numerical kernels are validated against an exact reference in
test_batched_ldsc_rg.py, so what this catches is the other failure mode -- a refactor that breaks
how the pieces are joined together.
"""

from pathlib import Path, PurePath
from typing import TypeVar

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.consolidate_ld_scores_task import (
    LD_SCORE_CHROM_COL,
    LD_SCORE_LD_SCORE_COL,
    LD_SCORE_M_5_50_COL,
    LD_SCORE_RSID_COL,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    BuildSlimProteinParquetTask,
    PppProteinFile,
)
from mecfs_bio.build_system.task.ppp_database.protein_sample_size_task import (
    PppProteinRef,
    PppProteinSampleSizeTask,
)
from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_GENE_CHROM_COL,
    ST3_GENE_END_COL,
    ST3_GENE_START_COL,
    ST3_OLINK_ID_COL,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import (
    ST3_VARIANT_SET_MISMATCH_ERR,
    PppProteinCrossTraitRgTask,
    PppRgConfig,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ppp_database_constants import (
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
    Oid,
    SynID,
)
from mecfs_bio.constants.ppp_ldsc_constants import (
    PPP_N_GENE_COL,
    PPP_N_OID_COL,
    PPP_N_SAMPLE_SIZE_COL,
    PPP_N_SYNID_COL,
    PPP_RG_GENE_COL,
    PPP_RG_H2_PROTEIN_COL,
    PPP_RG_H2_TRAIT_COL,
    PPP_RG_N_SNPS_COL,
    PPP_RG_OID_COL,
    PPP_RG_RG_COL,
    PPP_RG_RG_P_COL,
    PPP_RG_VARIANT_SET_COL,
    PPP_VARIANT_SET_ALL,
    PPP_VARIANT_SET_CIS_EXCLUDED,
    PppVariantSet,
)

# Small enough to build by hand, large enough for four jackknife blocks to be meaningful.
_N_SNPS = 400
_N_BLOCKS = 4
_TRAIT_N = 50_000
_PROTEIN_N = 30_000
# Simulated heritability, and the share of it the trait and the proteins hold in common. Real
# signal matters here: with pure noise every heritability lands at or below zero and rg is NaN
# by definition, which would leave the assertions unable to tell a working task from a broken
# one.
_H2 = 0.1
_SHARED_FRACTION = 0.6

# Two proteins sit on chr1, one on chr2. Only the chr1 pair appears in the gene-coordinate
# table, so the third has no cis window -- the case cis-excluded runs must skip.
_PROTEINS = (
    PppProteinFile(
        gene="AAA",
        uniprot="P00001",
        oid=Oid("OID00001"),
        version="v1",
        panel="Inflammation",
        synid=SynID("syn1"),
    ),
    PppProteinFile(
        gene="BBB",
        uniprot="P00002",
        oid=Oid("OID00002"),
        version="v1",
        panel="Oncology",
        synid=SynID("syn2"),
    ),
    PppProteinFile(
        gene="CCC",
        uniprot="P00003",
        oid=Oid("OID00003"),
        version="v1",
        panel="Neurology",
        synid=SynID("syn3"),
    ),
)
_PROTEINS_WITH_COORDS = (_PROTEINS[0], _PROTEINS[1])

_TaskT = TypeVar("_TaskT", bound=Task)


def _rsids() -> list[str]:
    return [f"rs{i}" for i in range(_N_SNPS)]


def _chroms() -> list[int]:
    # First half on chr1, second half on chr2; never chr6, so MHC exclusion is a no-op here.
    return [1] * (_N_SNPS // 2) + [2] * (_N_SNPS - _N_SNPS // 2)


def _positions() -> list[int]:
    return [1_000 + 1_000 * i for i in range(_N_SNPS)]


def _index_frame() -> pl.DataFrame:
    # A/G everywhere: unambiguous strand, so nothing is dropped as palindromic.
    return pl.DataFrame(
        {
            GWASLAB_CHROM_COL: _chroms(),
            GWASLAB_POS_COL: _positions(),
            GWASLAB_RSID_COL: _rsids(),
            GWASLAB_EFFECT_ALLELE_COL: ["A"] * _N_SNPS,
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G"] * _N_SNPS,
            PPP_INDEX_IS_STRAND_AMBIGUOUS_COL: [False] * _N_SNPS,
        }
    )


_LD_SCORES = np.random.default_rng(0).uniform(1.0, 10.0, _N_SNPS)
_M_5_50 = [5_000.0 if chrom == 1 else 4_000.0 for chrom in _chroms()]
_M_TOTAL = 9_000.0
# One shared per-variant genetic effect, which is what gives the trait and the proteins a
# genetic correlation to find.
_SHARED_EFFECT = np.random.default_rng(1).normal(0.0, 1.0, _N_SNPS)


def _ld_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            LD_SCORE_RSID_COL: _rsids(),
            LD_SCORE_CHROM_COL: _chroms(),
            LD_SCORE_LD_SCORE_COL: _LD_SCORES,
            LD_SCORE_M_5_50_COL: _M_5_50,
        }
    )


def _simulated_z(sample_size: int, seed: int) -> np.ndarray:
    """Z-scores under the LD-score model: variance grows with the LD score, as
    z ~ N(0, sqrt(1 + N*h2/M * l)), with part of the signal shared across traits."""
    rng = np.random.default_rng(seed)
    signal = np.sqrt(sample_size * _H2 / _M_TOTAL * _LD_SCORES)
    shared = np.sqrt(_SHARED_FRACTION) * _SHARED_EFFECT
    private = np.sqrt(1.0 - _SHARED_FRACTION) * rng.normal(0.0, 1.0, _N_SNPS)
    return signal * (shared + private) + rng.normal(0.0, 1.0, _N_SNPS)


def _protein_frame(seed: int) -> pl.DataFrame:
    """One protein's slim beta/se, in index row order."""
    se = np.random.default_rng(seed).uniform(0.01, 0.02, _N_SNPS)
    return pl.DataFrame(
        {
            GWASLAB_BETA_COL: _simulated_z(_PROTEIN_N, seed) * se,
            GWASLAB_SE_COL: se,
        }
    )


def _trait_frame() -> pl.DataFrame:
    """Trait sumstats covering every index variant, with the last one strand-flipped so the
    alignment step has something to negate."""
    se = np.random.default_rng(99).uniform(0.01, 0.02, _N_SNPS)
    effect = ["A"] * _N_SNPS
    other = ["G"] * _N_SNPS
    beta = _simulated_z(_TRAIT_N, 99) * se
    # The flipped variant reports the opposite allele, so its beta reverses sign too; alignment
    # has to undo exactly that.
    effect[-1], other[-1] = "G", "A"
    beta[-1] = -beta[-1]
    return pl.DataFrame(
        {
            GWASLAB_RSID_COL: _rsids(),
            GWASLAB_EFFECT_ALLELE_COL: effect,
            GWASLAB_NON_EFFECT_ALLELE_COL: other,
            GWASLAB_BETA_COL: beta,
            GWASLAB_SE_COL: se,
            GWASLAB_SAMPLE_SIZE_COLUMN: [float(_TRAIT_N)] * _N_SNPS,
        }
    )


def _sample_size_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            PPP_N_OID_COL: [str(protein.oid) for protein in _PROTEINS],
            PPP_N_GENE_COL: [protein.gene for protein in _PROTEINS],
            PPP_N_SYNID_COL: [str(protein.synid) for protein in _PROTEINS],
            PPP_N_SAMPLE_SIZE_COL: [_PROTEIN_N] * len(_PROTEINS),
        }
    )


def _gene_coords_frame() -> pl.DataFrame:
    """ST3 as extracted: gene windows for the first two proteins only."""
    return pl.DataFrame(
        {
            ST3_OLINK_ID_COL: [str(p.oid) for p in _PROTEINS_WITH_COORDS],
            ST3_GENE_CHROM_COL: ["1", "1"],
            ST3_GENE_START_COL: ["50000", "150000"],
            ST3_GENE_END_COL: ["60000", "160000"],
        }
    )


def _table_meta(asset_id: str) -> ResultTableMeta:
    return ResultTableMeta(
        id=AssetId(asset_id),
        trait="synthetic",
        project="test",
        sub_dir=PurePath("analysis"),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )


def _trait_meta(asset_id: str) -> GWASSummaryDataFileMeta:
    return GWASSummaryDataFileMeta(
        id=AssetId(asset_id),
        trait="synthetic_trait",
        project="synthetic_project",
        sub_dir="harmonized",
        project_path=PurePath(f"{asset_id}.parquet"),
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        extension=".parquet",
    )


def _build_task(
    tmp_path: Path, variant_set: PppVariantSet
) -> tuple[PppProteinCrossTraitRgTask, dict[AssetId, Asset]]:
    """The task under test, plus the fetch table mapping each dependency to a written parquet."""
    assets: dict[AssetId, Asset] = {}

    def register(task: _TaskT, frame: pl.DataFrame) -> _TaskT:
        path = tmp_path / f"{task.asset_id}.parquet"
        frame.write_parquet(path)
        assets[task.asset_id] = FileAsset(path)
        return task

    index_task = register(FakeTask(meta=_table_meta("index")), _index_frame())
    ld_scores_task = register(FakeTask(meta=_table_meta("ld_scores")), _ld_frame())
    trait_task = register(FakeTask(meta=_trait_meta("trait")), _trait_frame())
    gene_coords_task = register(
        FakeTask(meta=_table_meta("gene_coords")), _gene_coords_frame()
    )
    sample_size_task = PppProteinSampleSizeTask.create(
        asset_id="sample_sizes",
        protein_refs=tuple(
            PppProteinRef(oid=str(p.oid), gene=p.gene, synid=str(p.synid))
            for p in _PROTEINS
        ),
    )
    register(sample_size_task, _sample_size_frame())

    protein_tasks = tuple(
        register(
            BuildSlimProteinParquetTask.create(
                index_task=index_task,
                protein=protein,
                asset_id=f"slim_{protein.oid}",
                index_name="test",
            ),
            _protein_frame(seed=seed),
        )
        for seed, protein in enumerate(_PROTEINS)
    )

    task = PppProteinCrossTraitRgTask.create(
        asset_id="rg",
        trait_task=trait_task,
        protein_tasks=protein_tasks,
        index_task=index_task,
        ld_scores_task=ld_scores_task,
        sample_size_task=sample_size_task,
        gene_coords_task=(
            gene_coords_task if variant_set == PPP_VARIANT_SET_CIS_EXCLUDED else None
        ),
        config=PppRgConfig(
            variant_set=variant_set,
            n_blocks=_N_BLOCKS,
            batch_size=2,  # smaller than the protein count, so batching is exercised
            min_trait_snps=10,
        ),
    )
    return task, assets


def _run(tmp_path: Path, variant_set: PppVariantSet) -> pl.DataFrame:
    task, assets = _build_task(tmp_path, variant_set)
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir(exist_ok=True)
    result = task.execute(
        scratch_dir=scratch_dir, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_rg_task_all_variants(tmp_path: Path):
    table = _run(tmp_path, PPP_VARIANT_SET_ALL)

    # Every protein gets a row, identified from the protein task rather than any table lookup.
    assert table.height == len(_PROTEINS)
    assert set(table[PPP_RG_OID_COL]) == {str(p.oid) for p in _PROTEINS}
    assert set(table[PPP_RG_GENE_COL]) == {p.gene for p in _PROTEINS}
    assert table[PPP_RG_VARIANT_SET_COL].unique().to_list() == [PPP_VARIANT_SET_ALL]
    # Rows come out ordered by significance, so the strongest correlations read first.
    assert table[PPP_RG_RG_P_COL].to_list() == sorted(table[PPP_RG_RG_P_COL].to_list())

    # The trait heritability is estimated once and shared by every row.
    assert table[PPP_RG_H2_TRAIT_COL].n_unique() == 1
    # Every protein regressed on the full context, and the simulated shared signal comes back
    # as a real, positive genetic correlation.
    assert (table[PPP_RG_N_SNPS_COL] > 0).all()
    assert table[PPP_RG_RG_COL].is_finite().all()
    assert (table[PPP_RG_RG_COL] > 0).all()
    assert (table[PPP_RG_H2_PROTEIN_COL] > 0).all()
    assert (table[PPP_RG_H2_TRAIT_COL] > 0).all()


def test_rg_task_cis_excluded_skips_proteins_without_gene_coordinates(tmp_path: Path):
    table = _run(tmp_path, PPP_VARIANT_SET_CIS_EXCLUDED)

    # The third protein has no ST3 window, so there is nothing to exclude and no row.
    assert set(table[PPP_RG_OID_COL]) == {str(p.oid) for p in _PROTEINS_WITH_COORDS}
    assert table[PPP_RG_VARIANT_SET_COL].unique().to_list() == [
        PPP_VARIANT_SET_CIS_EXCLUDED
    ]
    # Dropping each gene's cis window leaves strictly fewer variants than the whole context.
    all_variants = _run(tmp_path, PPP_VARIANT_SET_ALL)
    context_snps = all_variants[PPP_RG_N_SNPS_COL].max()
    assert (table[PPP_RG_N_SNPS_COL] < context_snps).all()


def test_rg_task_requires_gene_coordinates_exactly_when_excluding_cis(tmp_path: Path):
    """The gene-coordinate dependency and the variant set have to agree: without it a
    cis-excluded run has no windows, and with it an all-variants run declares a dependency it
    never reads."""
    task, _ = _build_task(tmp_path, PPP_VARIANT_SET_CIS_EXCLUDED)
    with pytest.raises(AssertionError, match=ST3_VARIANT_SET_MISMATCH_ERR):
        PppProteinCrossTraitRgTask.create(
            asset_id="rg",
            trait_task=task.trait_task,
            protein_tasks=task.protein_tasks,
            index_task=task.index_task,
            ld_scores_task=task.ld_scores_task,
            sample_size_task=task.sample_size_task,
            gene_coords_task=None,
            config=PppRgConfig(variant_set=PPP_VARIANT_SET_CIS_EXCLUDED),
        )


def test_rg_task_deps_include_gene_coordinates_only_when_needed(tmp_path: Path):
    cis_task, _ = _build_task(tmp_path, PPP_VARIANT_SET_CIS_EXCLUDED)
    all_task, _ = _build_task(tmp_path, PPP_VARIANT_SET_ALL)
    assert cis_task.gene_coords_task in cis_task.deps
    assert all_task.gene_coords_task is None
    assert len(all_task.deps) == len(cis_task.deps) - 1
