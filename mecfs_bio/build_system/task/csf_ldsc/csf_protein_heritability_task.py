"""
Per-aptamer SNP heritability for the Western et al. 2024 CSF pQTL database via batched
LD-score regression.

This is the CSF analogue of PppProteinHeritabilityTask, deliberately reusing that
pipeline's shared machinery: every aptamer's slim file stores beta/se/N in the SAME CSF
variant-index row order, so the index<->LD-score alignment, LD scores, M and the
jackknife blocks are built once (build_batched_ldsc_context / BatchedLdscContext) and the
per-aptamer work is a cheap batched weighted regression (batched_h2).

Two differences from the PPP task:
  - Only all-variants heritability is produced (no cis-excluded set), so there is one
    row per aptamer and no gene-coordinate dependency.
  - CSF carries a per-variant N (unlike UKB-PPP's constant-per-protein N). The n_spread
    probe showed N is effectively constant over the HapMap3 regression SNPs (median->
    per-SNP h2 shift < ~0.03%), so we collapse it to one scalar per aptamer via the
    MEDIAN of the present-variant N (median, not the PPP equality assert, because CSF has
    a thin low-N tail).
"""

from __future__ import annotations

from pathlib import Path, PurePath

import numpy as np
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.batched_ldsc.batched_ldsc_context import (
    BatchedLdscContext,
    build_batched_ldsc_context,
)
from mecfs_bio.build_system.task.batched_ldsc.batched_ldsc_h2 import (
    DEFAULT_N_BLOCKS,
    batched_h2,
)
from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    BuildSlimCsfAptamerParquetTask,
)
from mecfs_bio.build_system.task.task_util import produces_dataframe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.csf_database_constants import (
    CSF_INDEX_IS_STRAND_AMBIGUOUS_COL,
)
from mecfs_bio.constants.csf_ldsc_constants import (
    CSF_H2_ANALYTE_COL,
    CSF_H2_GENE_SYMBOL_COL,
    CSF_H2_H2_COL,
    CSF_H2_H2_SE_COL,
    CSF_H2_INTERCEPT_COL,
    CSF_H2_LAMBDA_GC_COL,
    CSF_H2_MEAN_CHI2_COL,
    CSF_H2_N_BAR_COL,
    CSF_H2_N_SNPS_COL,
    CSF_H2_UNIPROT_COL,
    CSF_H2_VARIANT_SET_COL,
    CSF_VARIANT_SET_ALL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

logger = structlog.get_logger()

# Columns the shared context builder needs from the CSF variant index (row i aligns to
# row i of every slim file). is_strand_ambiguous is the same string PPP uses.
_INDEX_CONTEXT_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    CSF_INDEX_IS_STRAND_AMBIGUOUS_COL,
]

# Stable phrase in median_sample_size's failure message, shared with the test so the
# wording lives in one place.
NO_PRESENT_VARIANTS_ERR = "no present variants"


@frozen
class CsfHeritabilityConfig:
    drop_strand_ambiguous: bool = True
    exclude_mhc: bool = True
    n_blocks: int = DEFAULT_N_BLOCKS
    # Aptamers processed together per batched regression. Peak memory ~ n_snps * batch *
    # a few float64 arrays.
    batch_size: int = 50


def median_sample_size(n_at_context: np.ndarray, label: str) -> float:
    """The single N to use for an aptamer: the median of its per-variant N over the
    variants present at the context SNPs (absent variants are NaN and ignored). CSF N is
    effectively constant over the regression SNPs (see the n_spread probe), so the median
    is a faithful scalar; unlike the PPP equality assert it tolerates the thin low-N tail.
    label identifies the aptamer in the failure message."""
    finite = n_at_context[np.isfinite(n_at_context)]
    assert finite.size > 0, (
        f"aptamer {label} has {NO_PRESENT_VARIANTS_ERR} among its context SNPs; "
        "cannot recover a sample size"
    )
    return float(np.median(finite))


@frozen
class CsfProteinHeritabilityTask(GeneratingTask):
    """Compute all-variants LDSC heritability for every Western et al. 2024 CSF aptamer.

    aptamer_tasks: the slim per-aptamer beta/se/N tasks (all aligned to index_task).
    index_task: the shared ConstructCsfVariantIndexTask (variant identity / rsID / alleles).
    ld_scores_task: consolidated reference LD scores (ConsolidateLDScoresTask: CHR, SNP,
        L2, M_5_50).
    """

    meta: Meta
    aptamer_tasks: tuple[BuildSlimCsfAptamerParquetTask, ...]
    index_task: Task
    ld_scores_task: Task
    config: CsfHeritabilityConfig

    @property
    def deps(self) -> list[Task]:
        return [*self.aptamer_tasks, self.index_task, self.ld_scores_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        context = _build_context(
            self.index_task, self.ld_scores_task, self.config, fetch
        )
        logger.info(
            "built csf ldsc context",
            n_snps=context.n_snps,
            m=context.m,
            n_aptamers=len(self.aptamer_tasks),
        )

        rows: list[dict] = []
        tasks = list(self.aptamer_tasks)
        for start in range(0, len(tasks), self.config.batch_size):
            batch = tasks[start : start + self.config.batch_size]
            rows.extend(_process_batch(batch, context, self.config, fetch))

        table = pl.DataFrame(rows)
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet"
        table.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        aptamer_tasks: tuple[BuildSlimCsfAptamerParquetTask, ...],
        index_task: Task,
        ld_scores_task: Task,
        config: CsfHeritabilityConfig = CsfHeritabilityConfig(),
    ) -> CsfProteinHeritabilityTask:
        assert produces_dataframe(ld_scores_task), (
            f"ld_scores_task {ld_scores_task.asset_id} must produce a dataframe"
        )
        # Shared row order is the load-bearing invariant of the batched kernel: every
        # aptamer's slim file must be aligned to THIS index. It holds by construction
        # (BuildSlimCsfAptamerParquetTask aligns onto its index_task), but assert it so a
        # mismatched wiring fails at graph-build time rather than silently misaligning.
        for task in aptamer_tasks:
            assert task.index_task.asset_id == index_task.asset_id, (
                f"aptamer task {task.asset_id} is aligned to index "
                f"{task.index_task.asset_id}, not {index_task.asset_id}"
            )
        meta = ResultTableMeta(
            id=asset_id,
            trait="western_csf",
            project="csf_heritability",
            sub_dir=PurePath("analysis"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(
            meta=meta,
            aptamer_tasks=aptamer_tasks,
            index_task=index_task,
            ld_scores_task=ld_scores_task,
            config=config,
        )


def _build_context(
    index_task: Task,
    ld_scores_task: Task,
    config: CsfHeritabilityConfig,
    fetch: Fetch,
) -> BatchedLdscContext:
    index_df = (
        scan_dataframe_asset(
            fetch(index_task.asset_id),
            meta=index_task.meta,
            parquet_backend="polars",
        )
        .select(_INDEX_CONTEXT_COLUMNS)
        .collect()
        .to_polars()
    )
    return build_batched_ldsc_context(
        index_df,
        _read_table(ld_scores_task, fetch),
        drop_strand_ambiguous=config.drop_strand_ambiguous,
        exclude_mhc=config.exclude_mhc,
    )


def _read_table(task: Task, fetch: Fetch) -> pl.DataFrame:
    """A dependency's whole dataframe asset, in memory."""
    return (
        scan_dataframe_asset(
            fetch(task.asset_id),
            meta=task.meta,
            parquet_backend="polars",
        )
        .collect()
        .to_polars()
    )


def _process_batch(
    batch: list[BuildSlimCsfAptamerParquetTask],
    context: BatchedLdscContext,
    config: CsfHeritabilityConfig,
    fetch: Fetch,
) -> list[dict]:
    n_vec, chi2_cols = [], []
    for aptamer_task in batch:
        chi2, n = _read_aptamer_chi2(aptamer_task, context, fetch)
        chi2_cols.append(chi2)
        n_vec.append(n)
    chi2 = np.column_stack(chi2_cols)
    n_arr = np.asarray(n_vec)

    res = batched_h2(chi2, context.ld, n_arr, context.m, n_blocks=config.n_blocks)

    rows: list[dict] = []
    for j, aptamer_task in enumerate(batch):
        aptamer = aptamer_task.aptamer
        rows.append(
            {
                CSF_H2_ANALYTE_COL: aptamer.analyte,
                CSF_H2_UNIPROT_COL: aptamer.uniprot,
                CSF_H2_GENE_SYMBOL_COL: aptamer.entrez_gene_symbol,
                CSF_H2_VARIANT_SET_COL: CSF_VARIANT_SET_ALL,
                CSF_H2_H2_COL: float(res.h2[j]),
                CSF_H2_H2_SE_COL: float(res.h2_se[j]),
                CSF_H2_INTERCEPT_COL: float(res.intercept[j]),
                CSF_H2_MEAN_CHI2_COL: float(res.mean_chi2[j]),
                CSF_H2_LAMBDA_GC_COL: float(res.lambda_gc[j]),
                CSF_H2_N_SNPS_COL: int(res.n_snps[j]),
                CSF_H2_N_BAR_COL: float(n_arr[j]),
            }
        )
    return rows


def _read_aptamer_chi2(
    aptamer_task: BuildSlimCsfAptamerParquetTask,
    context: BatchedLdscContext,
    fetch: Fetch,
) -> tuple[np.ndarray, float]:
    """Return the aptamer's chi^2 at the context SNPs, plus its median sample size N over
    the present context variants."""
    frame = (
        scan_dataframe_asset(
            fetch(aptamer_task.asset_id),
            meta=aptamer_task.meta,
            parquet_backend="polars",
        )
        .select(GWASLAB_BETA_COL, GWASLAB_SE_COL, GWASLAB_SAMPLE_SIZE_COLUMN)
        .collect()
        .to_polars()
    )
    beta = frame[GWASLAB_BETA_COL].to_numpy().astype(float)[context.row_pos]
    se = frame[GWASLAB_SE_COL].to_numpy().astype(float)[context.row_pos]
    with np.errstate(invalid="ignore", divide="ignore"):
        chi2 = (beta / se) ** 2
    n_at_context = (
        frame[GWASLAB_SAMPLE_SIZE_COLUMN].to_numpy().astype(float)[context.row_pos]
    )
    n = median_sample_size(n_at_context, aptamer_task.asset_id)
    return chi2, n
