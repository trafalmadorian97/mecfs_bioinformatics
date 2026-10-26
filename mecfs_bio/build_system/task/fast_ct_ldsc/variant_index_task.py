"""
The shared regression-SNP index for fast all-pairs CT-LDSC.

The index is the HapMap3 SNP list (w_hm3.snplist: SNP, A1, A2) inner-joined on rsID to a
consolidated reference LD-score table (CHR, SNP, BP), with A1 as the effect allele and positions
taken from the LD scores (hg19 for the 1000 Genomes EUR reference). Rows are sorted by
(CHR, POS, rsID), and each row records whether the variant is strand-ambiguous.

Every downstream task derives the regression context from this index and the same LD scores with
build_ppp_ldsc_context(..., build="19"), which drops strand-ambiguous variants and the hg19
extended MHC. The context is a deterministic function of the two assets, so every aligned trait
and the all-pairs task see the same row set in the same order.
"""

from __future__ import annotations

from pathlib import Path, PurePath
from typing import Self

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.consolidate_ld_scores_task import (
    LD_SCORE_CHROM_COL,
    LD_SCORE_POS_COL,
    LD_SCORE_RSID_COL,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_rsid import (
    is_palindromic_expr,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_ldsc_context import (
    PppLdscContext,
    build_ppp_ldsc_context,
)
from mecfs_bio.build_system.task.task_util import produces_dataframe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.fast_ct_ldsc_constants import (
    HAPMAP3_SNPLIST_A1_COL,
    HAPMAP3_SNPLIST_A2_COL,
    HAPMAP3_SNPLIST_RSID_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.constants.ppp_database_constants import (
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
)

# Genome build of the index positions (taken from the 1000 Genomes EUR LD scores).
FAST_CT_LDSC_INDEX_BUILD: GenomeBuild = "19"

_CONTEXT_VARIANT_COLUMNS = [
    GWASLAB_RSID_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
]


@frozen(slots=True)
class FastCtLdscVariantIndexTask(Task):
    """Join the HapMap3 SNP list to the reference LD scores to form the regression-SNP index."""

    meta: Meta
    hapmap3_task: Task
    ld_scores_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.hapmap3_task, self.ld_scores_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        snplist = read_dataframe_task(self.hapmap3_task, fetch)
        ld_positions = read_dataframe_task(self.ld_scores_task, fetch).select(
            LD_SCORE_CHROM_COL, LD_SCORE_RSID_COL, LD_SCORE_POS_COL
        )
        index = _join_index(snplist, ld_positions)
        out_path = scratch_dir / "variant_index.parquet"
        index.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, asset_id: str, *, hapmap3_task: Task, ld_scores_task: Task
    ) -> Self:
        assert produces_dataframe(hapmap3_task), "hapmap3_task must produce a dataframe"
        assert produces_dataframe(ld_scores_task), (
            "ld_scores_task must produce a dataframe"
        )
        meta = ReferenceFileMeta(
            group="linkage_disequilibrium_score_reference_data",
            sub_group="fast_ct_ldsc",
            sub_folder=PurePath("processed"),
            id=AssetId(asset_id),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(meta=meta, hapmap3_task=hapmap3_task, ld_scores_task=ld_scores_task)


def _join_index(snplist: pl.DataFrame, ld_positions: pl.DataFrame) -> pl.DataFrame:
    return (
        snplist.select(
            pl.col(HAPMAP3_SNPLIST_RSID_COL).alias(GWASLAB_RSID_COL),
            pl.col(HAPMAP3_SNPLIST_A1_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
            pl.col(HAPMAP3_SNPLIST_A2_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        )
        .join(
            ld_positions.select(
                pl.col(LD_SCORE_CHROM_COL).cast(pl.Int64).alias(GWASLAB_CHROM_COL),
                pl.col(LD_SCORE_POS_COL).cast(pl.Int64).alias(GWASLAB_POS_COL),
                pl.col(LD_SCORE_RSID_COL).alias(GWASLAB_RSID_COL),
            ),
            on=GWASLAB_RSID_COL,
            how="inner",
        )
        .with_columns(
            is_palindromic_expr(
                GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
            ).alias(PPP_INDEX_IS_STRAND_AMBIGUOUS_COL)
        )
        .select(
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_RSID_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
            PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
        )
        .sort([GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_RSID_COL])
    )


def read_dataframe_task(task: Task, fetch: Fetch) -> pl.DataFrame:
    """A dependency's whole dataframe asset, in memory."""
    return (
        scan_dataframe_asset(
            fetch(task.asset_id), meta=task.meta, parquet_backend="polars"
        )
        .collect(backend="polars")
        .to_polars()
    )


@frozen(slots=True)
class IndexedLdscContext:
    """The regression context plus its variants (rsID, EA, NEA, CHR, POS) in context row order."""

    context: PppLdscContext
    variants: pl.DataFrame

    def __attrs_post_init__(self) -> None:
        assert self.variants.height == self.context.n_snps, (
            f"variants has {self.variants.height} rows, context has {self.context.n_snps}"
        )


def load_indexed_ldsc_context(
    index_task: FastCtLdscVariantIndexTask, fetch: Fetch
) -> IndexedLdscContext:
    index_df = read_dataframe_task(index_task, fetch)
    context = build_ppp_ldsc_context(
        index_df,
        read_dataframe_task(index_task.ld_scores_task, fetch),
        build=FAST_CT_LDSC_INDEX_BUILD,
    )
    variants = index_df.select(_CONTEXT_VARIANT_COLUMNS).gather(
        pl.Series(context.row_pos)
    )
    return IndexedLdscContext(context=context, variants=variants)
