"""Tests for CsfProteinHeritabilityTask.

Two levels: the median sample-size helper (CSF stores a per-variant N, collapsed to one
scalar per aptamer), and an end-to-end Task run over synthetic inputs that exercises the
wiring -- shared LDSC context, per-aptamer batching, N recovery and output assembly. The
numerical kernel (batched_h2) is validated against an exact reference in
test_batched_ldsc_h2.py, so this catches the other failure mode: a refactor that breaks
how the pieces are joined.
"""

from pathlib import Path, PurePath
from typing import TypeVar

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
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
from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    BuildSlimCsfAptamerParquetTask,
    CsfAptamerFile,
)
from mecfs_bio.build_system.task.csf_ldsc.csf_protein_heritability_task import (
    NO_PRESENT_VARIANTS_ERR,
    CsfHeritabilityConfig,
    CsfProteinHeritabilityTask,
    median_sample_size,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.csf_database_constants import (
    CSF_INDEX_IS_STRAND_AMBIGUOUS_COL,
    Analyte,
    GcstAccession,
    SeqId,
    UniProtId,
)
from mecfs_bio.constants.csf_ldsc_constants import (
    CSF_H2_ANALYTE_COL,
    CSF_H2_GENE_SYMBOL_COL,
    CSF_H2_H2_COL,
    CSF_H2_N_BAR_COL,
    CSF_H2_N_SNPS_COL,
    CSF_H2_UNIPROT_COL,
    CSF_H2_VARIANT_SET_COL,
    CSF_VARIANT_SET_ALL,
)
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


def test_median_sample_size_ignores_nan_absent_variants():
    n_at_context = np.array([np.nan, 3400.0, 3400.0, np.nan, 3400.0])
    assert median_sample_size(n_at_context, "APTAMER") == pytest.approx(3400.0)


def test_median_sample_size_tolerates_low_n_tail():
    # CSF's thin low-N tail must not trip the helper (unlike PPP's equality assert); the
    # median is robust to it.
    n_at_context = np.array([3160.0, 3400.0, 3400.0, 3400.0, 3400.0])
    assert median_sample_size(n_at_context, "APTAMER") == pytest.approx(3400.0)


def test_median_sample_size_raises_when_all_absent():
    with pytest.raises(AssertionError, match=NO_PRESENT_VARIANTS_ERR):
        median_sample_size(np.array([np.nan, np.nan]), "APTAMER")


# --------------------------------------------------------------------------------------
# End-to-end Task run over synthetic inputs.
# --------------------------------------------------------------------------------------
# Small enough to build by hand, large enough for the jackknife blocks to be meaningful.
_N_SNPS = 400
_N_BLOCKS = 4
_H2 = 0.1

_TaskT = TypeVar("_TaskT", bound=Task)

_APTAMERS = (
    CsfAptamerFile(
        analyte=Analyte("X1.1"),
        seq_id=SeqId("1-1"),
        accession=GcstAccession("GCST90400001"),
        uniprot=UniProtId("P00001"),
        entrez_gene_symbol="GENEA",
    ),
    CsfAptamerFile(
        analyte=Analyte("X2.2"),
        seq_id=SeqId("2-2"),
        accession=GcstAccession("GCST90400002"),
        uniprot=UniProtId("P00002"),
        entrez_gene_symbol="GENEB",
    ),
)
_APTAMER_N = (3400.0, 2500.0)


def _rsids() -> list[str]:
    return [f"rs{i}" for i in range(_N_SNPS)]


_LD_SCORES = np.random.default_rng(0).uniform(1.0, 10.0, _N_SNPS)
_M_TOTAL = 9_000.0


def _index_frame() -> pl.DataFrame:
    # A/G everywhere: unambiguous strand, so nothing is dropped as palindromic. All on
    # chr1, never chr6, so MHC exclusion is a no-op here.
    return pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1] * _N_SNPS,
            GWASLAB_POS_COL: [1_000 + 1_000 * i for i in range(_N_SNPS)],
            GWASLAB_EFFECT_ALLELE_COL: ["A"] * _N_SNPS,
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G"] * _N_SNPS,
            GWASLAB_RSID_COL: _rsids(),
            "EAF": np.random.default_rng(2).uniform(0.1, 0.9, _N_SNPS),
            CSF_INDEX_IS_STRAND_AMBIGUOUS_COL: [False] * _N_SNPS,
        }
    )


def _ld_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            LD_SCORE_RSID_COL: _rsids(),
            LD_SCORE_CHROM_COL: [1] * _N_SNPS,
            LD_SCORE_LD_SCORE_COL: _LD_SCORES,
            LD_SCORE_M_5_50_COL: [_M_TOTAL / _N_SNPS] * _N_SNPS,
        }
    )


def _aptamer_frame(sample_size: float, seed: int) -> pl.DataFrame:
    """One aptamer's slim beta/se/N in index row order, with a real heritability signal so
    h2 is estimable (pure noise lands at/below zero)."""
    rng = np.random.default_rng(seed)
    se = rng.uniform(0.01, 0.02, _N_SNPS)
    signal = np.sqrt(sample_size * _H2 / _M_TOTAL * _LD_SCORES)
    z = signal * rng.normal(0.0, 1.0, _N_SNPS) + rng.normal(0.0, 1.0, _N_SNPS)
    return pl.DataFrame(
        {
            GWASLAB_BETA_COL: (z * se).astype(np.float32),
            GWASLAB_SE_COL: se.astype(np.float32),
            GWASLAB_SAMPLE_SIZE_COLUMN: np.full(_N_SNPS, sample_size, dtype=np.float32),
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


def _build_task(
    tmp_path: Path,
) -> tuple[CsfProteinHeritabilityTask, dict[AssetId, Asset]]:
    assets: dict[AssetId, Asset] = {}

    def register(task: _TaskT, frame: pl.DataFrame) -> _TaskT:
        path = tmp_path / f"{task.asset_id}.parquet"
        frame.write_parquet(path)
        assets[task.asset_id] = FileAsset(path)
        return task

    index_task = register(FakeTask(meta=_table_meta("index")), _index_frame())
    ld_scores_task = register(FakeTask(meta=_table_meta("ld_scores")), _ld_frame())

    aptamer_tasks = tuple(
        register(
            BuildSlimCsfAptamerParquetTask.create(
                index_task=index_task,
                aptamer=aptamer,
                asset_id=f"slim_{aptamer.analyte}",
                index_name="test",
            ),
            _aptamer_frame(_APTAMER_N[i], seed=i + 1),
        )
        for i, aptamer in enumerate(_APTAMERS)
    )

    task = CsfProteinHeritabilityTask.create(
        asset_id="csf_h2",
        aptamer_tasks=aptamer_tasks,
        index_task=index_task,
        ld_scores_task=ld_scores_task,
        config=CsfHeritabilityConfig(
            n_blocks=_N_BLOCKS,
            batch_size=1,  # smaller than the aptamer count, so batching is exercised
        ),
    )
    return task, assets


def test_task_produces_one_row_per_aptamer(tmp_path: Path):
    task, assets = _build_task(tmp_path)
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir(exist_ok=True)
    result = task.execute(
        scratch_dir=scratch_dir, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, FileAsset)
    table = pl.read_parquet(result.path)

    assert table.height == len(_APTAMERS)
    assert table[CSF_H2_ANALYTE_COL].to_list() == ["X1.1", "X2.2"]
    assert table[CSF_H2_UNIPROT_COL].to_list() == ["P00001", "P00002"]
    assert table[CSF_H2_GENE_SYMBOL_COL].to_list() == ["GENEA", "GENEB"]
    # N is recovered as the median per-variant N (constant here).
    assert table[CSF_H2_N_BAR_COL].to_list() == list(_APTAMER_N)
    # Forward-compat variant-set column is present and constant.
    assert table[CSF_H2_VARIANT_SET_COL].unique().to_list() == [CSF_VARIANT_SET_ALL]
    assert table[CSF_H2_H2_COL].null_count() == 0
    # Every SNP survives the (no-op) filters here, so all are used in the regression.
    assert table[CSF_H2_N_SNPS_COL].to_list() == [_N_SNPS, _N_SNPS]


def test_create_rejects_aptamer_aligned_to_other_index(tmp_path: Path):
    _task, _assets = _build_task(tmp_path)
    other_index = FakeTask(meta=_table_meta("other_index"))
    ld_scores_task = FakeTask(meta=_table_meta("ld_scores"))
    aptamer_task = BuildSlimCsfAptamerParquetTask.create(
        index_task=other_index,
        aptamer=_APTAMERS[0],
        asset_id="slim_mismatched",
        index_name="test",
    )
    # The aptamer is aligned to other_index but we wire a different index_task.
    with pytest.raises(AssertionError):
        CsfProteinHeritabilityTask.create(
            asset_id="csf_h2",
            aptamer_tasks=(aptamer_task,),
            index_task=FakeTask(meta=_table_meta("index")),
            ld_scores_task=ld_scores_task,
        )
