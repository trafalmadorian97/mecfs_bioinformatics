"""
Per-protein SNP heritability for the deCODE SomaScan database via batched LD-score regression,
alongside the UKB-PPP task and sharing its regression context and kernel.

For every protein this produces two heritabilities: one on all regression variants, and one
excluding cis variants (within a window of the protein's gene, from the MAGMA GRCh38 table via
DecodeGeneCoordsReader). Proteins without gene coordinates get only the all-variants row.

deCODE's N varies by variant (imputation and per-sample missingness), so each protein's N is the
median of its per-variant N over the regression variants (n_bar), with the 5th and 95th
percentiles (n_p05, n_p95) reported to show the spread.
"""

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
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.decode_database.build_decode_somascan_protein_task import (
    DECODE_SLIM_MEMBER,
    BuildDecodeSomascanProteinTask,
)
from mecfs_bio.build_system.task.directory_member_dataframe import (
    scan_directory_member,
)
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_h2 import (
    BatchedH2Result,
    batched_h2,
)
from mecfs_bio.build_system.task.ppp_ldsc.heritability_common import (
    PppHeritabilityConfig,
    build_cis_exclude,
    build_context,
    read_table,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_ldsc_context import PppLdscContext
from mecfs_bio.build_system.task.protein_coloc.decode_gene_coords_reader import (
    DecodeGeneCoordsReader,
)
from mecfs_bio.build_system.task.task_util import produces_dataframe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_H2_N_P05_COL,
    DECODE_H2_N_P95_COL,
    DECODE_H2_PROTEIN_ID_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ppp_ldsc_constants import (
    PPP_H2_GENE_COL,
    PPP_H2_H2_COL,
    PPP_H2_H2_SE_COL,
    PPP_H2_INTERCEPT_COL,
    PPP_H2_LAMBDA_GC_COL,
    PPP_H2_MEAN_CHI2_COL,
    PPP_H2_N_BAR_COL,
    PPP_H2_N_SNPS_COL,
    PPP_H2_VARIANT_SET_COL,
    PPP_VARIANT_SET_ALL,
    PPP_VARIANT_SET_CIS_EXCLUDED,
    PppVariantSet,
)
from mecfs_bio.constants.protein_coloc_constants import ProteinId
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

logger = structlog.get_logger()


@frozen(slots=True)
class _ProteinN:
    """Median, 5th and 95th percentile of a protein's per-variant N over the context."""

    median: float
    p05: float
    p95: float


@frozen(slots=True)
class _ProteinStats:
    chi2: np.ndarray
    n: _ProteinN


@frozen(slots=True)
class DecodeProteinHeritabilityTask(GeneratingTask):
    """
    protein_tasks: the per-protein deCODE directories, whose slim files are aligned to
        index_task; each task's protein gives the seq id, gene symbol and Entrez ids.
    index_task: the group's ConstructDecodeVariantIndexTask.
    ld_scores_task: consolidated reference LD scores (ConsolidateLDScoresTask).
    gene_coords_task: the MAGMA GRCh38 gene-location table.
    """

    meta: Meta
    protein_tasks: tuple[BuildDecodeSomascanProteinTask, ...]
    index_task: Task
    ld_scores_task: Task
    gene_coords_task: Task
    config: PppHeritabilityConfig

    @property
    def deps(self) -> list[Task]:
        return [
            *self.protein_tasks,
            self.index_task,
            self.ld_scores_task,
            self.gene_coords_task,
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        context = build_context(
            self.index_task, self.ld_scores_task, self.config, fetch
        )
        gene_coords = DecodeGeneCoordsReader(
            proteins=tuple(task.protein for task in self.protein_tasks)
        )(read_table(self.gene_coords_task, fetch))
        logger.info(
            "built decode ldsc context",
            n_snps=context.n_snps,
            m=context.m,
            n_proteins=len(self.protein_tasks),
        )
        tasks = list(self.protein_tasks)
        batches = [
            tasks[start : start + self.config.batch_size]
            for start in range(0, len(tasks), self.config.batch_size)
        ]
        table = pl.DataFrame(
            [
                row
                for batch in batches
                for row in _process_batch(
                    batch, context, gene_coords, self.config, fetch
                )
            ]
        )
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet"
        table.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        protein_tasks: tuple[BuildDecodeSomascanProteinTask, ...],
        index_task: Task,
        ld_scores_task: Task,
        gene_coords_task: Task,
        config: PppHeritabilityConfig = PppHeritabilityConfig(),
    ) -> "DecodeProteinHeritabilityTask":
        assert produces_dataframe(ld_scores_task), (
            f"ld_scores_task {ld_scores_task.asset_id} must produce a dataframe"
        )
        assert produces_dataframe(gene_coords_task), (
            f"gene_coords_task {gene_coords_task.asset_id} must produce a dataframe"
        )
        traits = {task.meta.trait for task in protein_tasks}
        assert len(traits) == 1, f"protein tasks span several groups: {traits}"
        return cls(
            meta=ResultTableMeta(
                id=asset_id,
                trait=traits.pop(),
                project="decode_heritability",
                sub_dir=PurePath("analysis"),
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            protein_tasks=protein_tasks,
            index_task=index_task,
            ld_scores_task=ld_scores_task,
            gene_coords_task=gene_coords_task,
            config=config,
        )


def _read_protein(
    task: BuildDecodeSomascanProteinTask, context: PppLdscContext, fetch: Fetch
) -> _ProteinStats:
    """The protein's chi^2 at the context variants and the spread of its N there."""
    frame = (
        scan_directory_member(fetch(task.asset_id), DECODE_SLIM_MEMBER)
        .select(GWASLAB_BETA_COL, GWASLAB_SE_COL, GWASLAB_SAMPLE_SIZE_COLUMN)
        .collect()
        .to_polars()
    )
    beta = frame[GWASLAB_BETA_COL].to_numpy().astype(float)[context.row_pos]
    se = frame[GWASLAB_SE_COL].to_numpy().astype(float)[context.row_pos]
    with np.errstate(invalid="ignore", divide="ignore"):
        chi2 = (beta / se) ** 2
    n = frame[GWASLAB_SAMPLE_SIZE_COLUMN].to_numpy().astype(float)[context.row_pos]
    finite = n[np.isfinite(n)]
    assert finite.size > 0, f"protein {task.asset_id} has no N at the context variants"
    return _ProteinStats(
        chi2=chi2,
        n=_ProteinN(
            median=float(np.median(finite)),
            p05=float(np.quantile(finite, 0.05)),
            p95=float(np.quantile(finite, 0.95)),
        ),
    )


def _process_batch(
    batch: list[BuildDecodeSomascanProteinTask],
    context: PppLdscContext,
    gene_coords: dict[ProteinId, GenomicInterval],
    config: PppHeritabilityConfig,
    fetch: Fetch,
) -> list[dict]:
    stats = [_read_protein(task, context, fetch) for task in batch]
    protein_ids = [ProteinId(task.protein.seq_id) for task in batch]
    chi2 = np.column_stack([s.chi2 for s in stats])
    n_arr = np.asarray([s.n.median for s in stats])
    all_res = batched_h2(chi2, context.ld, n_arr, context.m, n_blocks=config.n_blocks)
    cis_res = batched_h2(
        chi2,
        context.ld,
        n_arr,
        context.m,
        n_blocks=config.n_blocks,
        exclude=build_cis_exclude(
            protein_ids, context, gene_coords, config.cis_window_bp
        ),
    )
    rows: list[dict] = []
    for j, task in enumerate(batch):
        rows.append(_make_row(task, stats[j].n, PPP_VARIANT_SET_ALL, all_res, j))
        if protein_ids[j] in gene_coords:
            rows.append(
                _make_row(task, stats[j].n, PPP_VARIANT_SET_CIS_EXCLUDED, cis_res, j)
            )
    return rows


def _make_row(
    task: BuildDecodeSomascanProteinTask,
    n: _ProteinN,
    variant_set: PppVariantSet,
    res: BatchedH2Result,
    j: int,
) -> dict:
    return {
        DECODE_H2_PROTEIN_ID_COL: str(task.protein.seq_id),
        PPP_H2_GENE_COL: task.protein.gene_symbol,
        PPP_H2_VARIANT_SET_COL: variant_set,
        PPP_H2_H2_COL: float(res.h2[j]),
        PPP_H2_H2_SE_COL: float(res.h2_se[j]),
        PPP_H2_INTERCEPT_COL: float(res.intercept[j]),
        PPP_H2_MEAN_CHI2_COL: float(res.mean_chi2[j]),
        PPP_H2_LAMBDA_GC_COL: float(res.lambda_gc[j]),
        PPP_H2_N_SNPS_COL: int(res.n_snps[j]),
        PPP_H2_N_BAR_COL: n.median,
        DECODE_H2_N_P05_COL: n.p05,
        DECODE_H2_N_P95_COL: n.p95,
    }
