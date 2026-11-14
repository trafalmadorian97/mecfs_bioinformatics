"""
Docs-facing table of the colocalizations of one trait with proteins.

Keeps the tested (locus, protein) pairs of a TraitProteinColocTask result whose PP.H4 exceeds
h4_threshold, sorted by PP.H4 descending. Column names are short so the table fits a docs page:
gene, protein (aptamer/protein id), locus (the locus id, chr:start-end in GRCh38), lead (the
trait lead as chr:pos:NEA:EA), cis (cis/trans status), pp_h4, q, fdr_sig, and:
- trait_p: p-value of the trait lead variant of the locus (over the whole trait window).
- prot_p: p-value of the protein lead variant within the aligned window.
- same_sign: whether the trait and protein effects have the same sign at the protein lead
  variant (the aligned variant with the largest |BETA_PROTEIN / SE_PROTEIN|), read from the
  pair's saved window, where both betas are relative to the same effect allele.
- n_snp: the number of aligned variants passed to coloc.abf. The prior odds of H3 against H4
  grow roughly linearly with it.

q is the Bayesian false discovery rate (Newton et al. 2004) of calling colocalization for
this pair and every tested pair ranked above it. With the tested pairs sorted by PP.H4
descending, the k-th pair gets

    q_k = (1/k) * sum_{i <= k} (1 - PP.H4_i),

the posterior expected fraction of false colocalizations among the top k. This needs no
independence between pairs, since expectations add regardless. Because 1 - PP.H4 does not
decrease down the ranking, q_k does not decrease either, so fdr_sig (q_k <= fdr_alpha)
marks exactly the largest top-k set whose expected FDR is at most fdr_alpha. The q-values are
only as well calibrated as PP.H4 itself, which depends on the coloc priors.
Written as zstd parquet with byte-stream-split floats and Int32 integers, for the data_table
docs macro. An empty table is a valid output.
"""

from pathlib import Path, PurePath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    float_column_names,
    write_parquet_table,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    CIS_STATUS_COL,
    DISPLAY_CIS_STATUS_COL,
    DISPLAY_FDR_SIGNIFICANT_COL,
    DISPLAY_GENE_COL,
    DISPLAY_LEAD_VARIANT_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_N_SNP_COL,
    DISPLAY_PP_H4_COL,
    DISPLAY_PROTEIN_COL,
    DISPLAY_PROTEIN_P_COL,
    DISPLAY_Q_VALUE_COL,
    DISPLAY_SAME_SIGN_COL,
    DISPLAY_TRAIT_P_COL,
    GENE_COL,
    LOCUS_CHR_COL,
    LOCUS_ID_COL,
    LOCUS_LEAD_EA_COL,
    LOCUS_LEAD_MLOG10P_COL,
    LOCUS_LEAD_NEA_COL,
    LOCUS_LEAD_POS_COL,
    N_ALIGNED_COL,
    PP_H4_COL,
    PROTEIN_ID_COL,
    PROTEIN_LEAD_MLOG10P_COL,
    RAW_RESULTS_FILENAME,
    SE_PROTEIN_COL,
    STATUS_COL,
    STATUS_TESTED,
    WINDOWS_DIRNAME,
    pair_file_stem,
)

DISPLAY_SCHEMA = pl.Schema(
    {
        DISPLAY_GENE_COL: pl.String,
        DISPLAY_PROTEIN_COL: pl.String,
        DISPLAY_LOCUS_COL: pl.String,
        DISPLAY_LEAD_VARIANT_COL: pl.String,
        DISPLAY_CIS_STATUS_COL: pl.String,
        DISPLAY_PP_H4_COL: pl.Float64,
        DISPLAY_Q_VALUE_COL: pl.Float64,
        DISPLAY_FDR_SIGNIFICANT_COL: pl.Boolean,
        DISPLAY_TRAIT_P_COL: pl.Float64,
        DISPLAY_PROTEIN_P_COL: pl.Float64,
        DISPLAY_SAME_SIGN_COL: pl.Boolean,
        DISPLAY_N_SNP_COL: pl.Int32,
    }
)


def coloc_result_dir(fetch: Fetch, coloc_task: Task) -> Path:
    asset = fetch(coloc_task.asset_id)
    assert isinstance(asset, DirectoryAsset), (
        f"expected {coloc_task.asset_id} to be a DirectoryAsset"
    )
    return asset.path


def read_raw_results(fetch: Fetch, coloc_task: Task) -> pl.DataFrame:
    return pl.read_parquet(coloc_result_dir(fetch, coloc_task) / RAW_RESULTS_FILENAME)


def _lead_signs_agree(windows_dir: Path, locus_id: str, protein_id: str) -> bool:
    """Whether trait and protein betas share a sign at the protein lead of the pair's window."""
    window = pl.read_parquet(
        windows_dir / f"{pair_file_stem(locus_id, protein_id)}.parquet"
    )
    lead_index = (window[BETA_PROTEIN_COL] / window[SE_PROTEIN_COL]).abs().arg_max()
    assert lead_index is not None, f"empty window for {locus_id} x {protein_id}"
    lead = window.row(lead_index, named=True)
    return (lead[BETA_TRAIT_COL] > 0) == (lead[BETA_PROTEIN_COL] > 0)


def _same_sign_column(rows: pl.DataFrame, windows_dir: Path) -> pl.Series:
    return pl.Series(
        DISPLAY_SAME_SIGN_COL,
        [
            _lead_signs_agree(windows_dir, locus_id, protein_id)
            for locus_id, protein_id in rows.select(
                LOCUS_ID_COL, PROTEIN_ID_COL
            ).iter_rows()
        ],
        dtype=pl.Boolean,
    )


def ranked_display_loci(display: pl.DataFrame, max_loci: int) -> list[str]:
    """Distinct display loci ordered by their best PP.H4 (descending), at most max_loci."""
    return (
        display.group_by(DISPLAY_LOCUS_COL)
        .agg(pl.col(DISPLAY_PP_H4_COL).max())
        .sort([DISPLAY_PP_H4_COL, DISPLAY_LOCUS_COL], descending=[True, False])
        .head(max_loci)[DISPLAY_LOCUS_COL]
        .to_list()
    )


def display_rows(
    raw: pl.DataFrame, windows_dir: Path, h4_threshold: float, fdr_alpha: float
) -> pl.DataFrame:
    lead_variant = pl.format(
        "chr{}:{}:{}:{}",
        LOCUS_CHR_COL,
        LOCUS_LEAD_POS_COL,
        LOCUS_LEAD_NEA_COL,
        LOCUS_LEAD_EA_COL,
    )
    q_value = (1 - pl.col(PP_H4_COL)).cum_sum() / pl.int_range(1, pl.len() + 1)
    hits = (
        raw.filter(pl.col(STATUS_COL) == STATUS_TESTED)
        .sort([PP_H4_COL, PROTEIN_ID_COL], descending=[True, False])
        .with_columns(q_value.alias(DISPLAY_Q_VALUE_COL))
        .filter(pl.col(PP_H4_COL) > h4_threshold)
    )
    return (
        hits.with_columns(_same_sign_column(hits, windows_dir))
        .select(
            pl.col(GENE_COL).alias(DISPLAY_GENE_COL),
            pl.col(PROTEIN_ID_COL).alias(DISPLAY_PROTEIN_COL),
            pl.col(LOCUS_ID_COL).alias(DISPLAY_LOCUS_COL),
            lead_variant.alias(DISPLAY_LEAD_VARIANT_COL),
            pl.col(CIS_STATUS_COL).alias(DISPLAY_CIS_STATUS_COL),
            pl.col(PP_H4_COL).alias(DISPLAY_PP_H4_COL),
            DISPLAY_Q_VALUE_COL,
            (pl.col(DISPLAY_Q_VALUE_COL) <= fdr_alpha).alias(
                DISPLAY_FDR_SIGNIFICANT_COL
            ),
            (10 ** -pl.col(LOCUS_LEAD_MLOG10P_COL)).alias(DISPLAY_TRAIT_P_COL),
            (10 ** -pl.col(PROTEIN_LEAD_MLOG10P_COL)).alias(DISPLAY_PROTEIN_P_COL),
            DISPLAY_SAME_SIGN_COL,
            pl.col(N_ALIGNED_COL).alias(DISPLAY_N_SNP_COL),
        )
        .cast(DISPLAY_SCHEMA)
    )


@frozen(slots=True)
class ColocDisplayTableTask(Task):
    meta: Meta
    coloc_task: Task
    h4_threshold: float
    fdr_alpha: float

    @property
    def deps(self) -> list[Task]:
        return [self.coloc_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        coloc_dir = coloc_result_dir(fetch, self.coloc_task)
        table = display_rows(
            pl.read_parquet(coloc_dir / RAW_RESULTS_FILENAME),
            windows_dir=coloc_dir / WINDOWS_DIRNAME,
            h4_threshold=self.h4_threshold,
            fdr_alpha=self.fdr_alpha,
        )
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet"
        arrow_table = table.to_arrow()
        write_parquet_table(
            table=arrow_table,
            out_path=out_path,
            compression="zstd",
            compression_level=None,
            byte_stream_split_columns=float_column_names(arrow_table),
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        coloc_task: Task,
        h4_threshold: float = 0.8,
        fdr_alpha: float = 0.05,
    ) -> "ColocDisplayTableTask":
        """coloc_task: a TraitProteinColocTask (any Task producing its result directory)."""
        assert 0 < h4_threshold < 1, f"h4_threshold must be in (0, 1): {h4_threshold}"
        assert 0 < fdr_alpha < 1, f"fdr_alpha must be in (0, 1): {fdr_alpha}"
        source_meta = coloc_task.meta
        assert isinstance(source_meta, ResultDirectoryMeta), (
            f"coloc_task must carry ResultDirectoryMeta, got {type(source_meta).__name__}"
        )
        return cls(
            meta=ResultTableMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=PurePath("analysis"),
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            coloc_task=coloc_task,
            h4_threshold=h4_threshold,
            fdr_alpha=fdr_alpha,
        )
