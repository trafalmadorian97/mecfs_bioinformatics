"""
Table of genome-wide significant hit segments for every UKB-PPP protein.

This screening table tells the colocalization Task which proteins have a pQTL near a trait
locus, without touching the dense per-protein summary statistics. It is computed from the slim
per-protein beta/se files over a shared variant index (for example the HapMap3 index), one
protein at a time.

For each protein and chromosome, the variants with p < p_threshold are sorted by position, and a
new segment starts whenever the gap to the previous significant variant exceeds merge_gap_bp.
This is a distance merge, not LD clumping: SEG_START and SEG_END are the positions (GRCh38) of
the first and last significant variant of the segment, so overlap with a trait window is exact.
Each segment reports its lead (min-p) variant and the number of significant variants in it.

p is the two-sided normal p-value of z = BETA / SE, kept in log space
(MLOG10P = -(log 2 + log sf(|z|)) / log 10) so that very strong pQTLs do not underflow.
"""

import math
from pathlib import Path, PurePath

import numpy as np
import polars as pl
from attrs import frozen
from scipy.stats import norm

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    BuildSlimProteinParquetTask,
)
from mecfs_bio.build_system.task.task_util import produces_dataframe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_MLOG10P_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    HIT_CHR_COL,
    HIT_GENE_COL,
    HIT_LEAD_MLOG10P_COL,
    HIT_LEAD_POS_COL,
    HIT_N_SIG_VARIANTS_COL,
    HIT_PROTEIN_ID_COL,
    HIT_SEG_END_COL,
    HIT_SEG_START_COL,
)

_SEGMENT_COL = "_segment"

HIT_TABLE_SCHEMA = pl.Schema(
    {
        HIT_PROTEIN_ID_COL: pl.String,
        HIT_GENE_COL: pl.String,
        HIT_CHR_COL: pl.Int32,
        HIT_SEG_START_COL: pl.Int64,
        HIT_SEG_END_COL: pl.Int64,
        HIT_LEAD_POS_COL: pl.Int64,
        HIT_LEAD_MLOG10P_COL: pl.Float64,
        HIT_N_SIG_VARIANTS_COL: pl.Int32,
    }
)


def two_sided_mlog10p(z: np.ndarray) -> np.ndarray:
    """-log10 of the two-sided normal p-value of z, computed in log space."""
    return -(math.log(2) + norm.logsf(np.abs(z))) / math.log(10)


def _read_index_positions(index_task: Task, fetch: Fetch) -> pl.DataFrame:
    return (
        scan_dataframe_asset(
            fetch(index_task.asset_id), meta=index_task.meta, parquet_backend="polars"
        )
        .to_native()
        .select(
            pl.col(GWASLAB_CHROM_COL).cast(pl.Int32),
            pl.col(GWASLAB_POS_COL).cast(pl.Int64),
        )
        .collect()
    )


def _significant_variants(
    protein_task: BuildSlimProteinParquetTask,
    index: pl.DataFrame,
    min_mlog10p: float,
    fetch: Fetch,
) -> pl.DataFrame:
    """Index CHR/POS plus MLOG10P of the protein's variants passing the threshold."""
    slim = (
        scan_dataframe_asset(
            fetch(protein_task.asset_id),
            meta=protein_task.meta,
            parquet_backend="polars",
        )
        .to_native()
        .select(GWASLAB_BETA_COL, GWASLAB_SE_COL)
        .collect()
    )
    assert slim.height == index.height, (
        f"{protein_task.asset_id} has {slim.height} rows but the index has {index.height}"
    )
    beta = slim[GWASLAB_BETA_COL].to_numpy().astype(np.float64)
    se = slim[GWASLAB_SE_COL].to_numpy().astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        mlog10p = two_sided_mlog10p(beta / se)
    keep = np.isfinite(mlog10p) & (mlog10p > min_mlog10p)
    return index.filter(pl.Series(keep)).with_columns(
        pl.Series(GWASLAB_MLOG10P_COL, mlog10p[keep])
    )


def merge_into_segments(significant: pl.DataFrame, merge_gap_bp: int) -> pl.DataFrame:
    """Distance-merge one protein's significant variants (CHR, POS, MLOG10P) into segments."""
    pos = pl.col(GWASLAB_POS_COL)
    return (
        significant.sort(GWASLAB_CHROM_COL, GWASLAB_POS_COL)
        .with_columns(
            (pos.diff().over(GWASLAB_CHROM_COL).fill_null(0) > merge_gap_bp)
            .cum_sum()
            .over(GWASLAB_CHROM_COL)
            .alias(_SEGMENT_COL)
        )
        .group_by(GWASLAB_CHROM_COL, _SEGMENT_COL)
        .agg(
            pos.min().alias(HIT_SEG_START_COL),
            pos.max().alias(HIT_SEG_END_COL),
            pos.sort_by(pl.col(GWASLAB_MLOG10P_COL), descending=True)
            .first()
            .alias(HIT_LEAD_POS_COL),
            pl.col(GWASLAB_MLOG10P_COL).max().alias(HIT_LEAD_MLOG10P_COL),
            pl.len().alias(HIT_N_SIG_VARIANTS_COL),
        )
        .rename({GWASLAB_CHROM_COL: HIT_CHR_COL})
        .drop(_SEGMENT_COL)
    )


def _protein_segments(
    protein_task: BuildSlimProteinParquetTask,
    index: pl.DataFrame,
    min_mlog10p: float,
    merge_gap_bp: int,
    fetch: Fetch,
) -> pl.DataFrame:
    segments = merge_into_segments(
        _significant_variants(protein_task, index, min_mlog10p, fetch), merge_gap_bp
    )
    return segments.with_columns(
        pl.lit(str(protein_task.protein.oid)).alias(HIT_PROTEIN_ID_COL),
        pl.lit(protein_task.protein.gene).alias(HIT_GENE_COL),
    )


@frozen(slots=True)
class PppProteinHitTableTask(GeneratingTask):
    """
    protein_tasks: the slim per-protein beta/se tasks, aligned to index_task.
    index_task: the shared variant index (CHR, POS in GRCh38).
    p_threshold: a variant is significant when p < p_threshold.
    merge_gap_bp: significant variants closer than this join one segment.
    """

    meta: Meta
    protein_tasks: tuple[BuildSlimProteinParquetTask, ...]
    index_task: Task
    p_threshold: float
    merge_gap_bp: int

    @property
    def deps(self) -> list[Task]:
        return [self.index_task, *self.protein_tasks]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        index = _read_index_positions(self.index_task, fetch)
        min_mlog10p = -math.log10(self.p_threshold)
        segments = [
            _protein_segments(
                protein_task, index, min_mlog10p, self.merge_gap_bp, fetch
            ).select(HIT_TABLE_SCHEMA.names())
            for protein_task in self.protein_tasks
        ]
        table = (
            pl.concat(segments)
            .cast(HIT_TABLE_SCHEMA)
            .sort(HIT_PROTEIN_ID_COL, HIT_CHR_COL, HIT_SEG_START_COL)
            if segments
            else HIT_TABLE_SCHEMA.to_frame()
        )
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet"
        table.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        protein_tasks: tuple[BuildSlimProteinParquetTask, ...],
        index_task: Task,
        p_threshold: float = 5e-8,
        merge_gap_bp: int = 250_000,
    ) -> "PppProteinHitTableTask":
        assert 0 < p_threshold < 1, f"p_threshold must be in (0, 1): {p_threshold}"
        assert merge_gap_bp >= 0, f"merge_gap_bp must be non-negative: {merge_gap_bp}"
        assert produces_dataframe(index_task), (
            f"index_task {index_task.asset_id} must produce a dataframe"
        )
        meta = ResultTableMeta(
            id=AssetId(asset_id),
            trait="ukbb_ppp",
            project="ppp_hit_segments",
            sub_dir=PurePath("analysis"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(
            meta=meta,
            protein_tasks=protein_tasks,
            index_task=index_task,
            p_threshold=p_threshold,
            merge_gap_bp=merge_gap_bp,
        )
