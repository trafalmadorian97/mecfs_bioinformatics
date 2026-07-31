"""
Normalized common-1kg membership list for the PPP variant index (common-1kg mode).

Takes the filtered common biallelic SNVs (FilterCommon1kgVariantsTask) and attaches
rsIDs by an allele-aware join to the gwaslab 1000G-dbSNP151 hg38 table, then emits
the membership contract (CHR, POS[hg38], EA, NEA, rsID). rsID is nullable: variants
absent from the dbSNP table keep a null rsID. Allele orientation is arbitrary (the
index adopts the template protein's orientation); EA=ALT, NEA=REF here.
"""

from pathlib import Path, PurePath

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
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_util import (
    gwaslab_download_ref_if_missing,
)
from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key
from mecfs_bio.build_system.task.ppp_database.byte_stream_split_parquet import (
    write_byte_stream_split_parquet,
)
from mecfs_bio.build_system.task.ppp_database.filter_common_1kg_variants_task import (
    ONEKG_ALT_COL,
    ONEKG_REF_COL,
)
from mecfs_bio.build_system.task.ppp_database.hapmap3_membership_task import (
    MEMBERSHIP_COLUMNS,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_1KG_DBSNP151_HG38_AUTO_NAME,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.constants.ppp_database_constants import PPP_INDEX_ALLELE_KEY_COL

_DBSNP_SEPARATOR = "\t"


def attach_rsid(filtered: pl.DataFrame, dbsnp: pl.DataFrame) -> pl.DataFrame:
    """Left-join rsIDs from a dbSNP table (columns CHR, POS, EA, NEA, rsID) onto the
    filtered 1kg variants (columns CHR, POS, REF, ALT, AF) via an allele-aware key,
    returning the membership contract (CHR, POS, EA=ALT, NEA=REF, rsID)."""
    filtered_keyed = filtered.with_columns(
        unordered_allele_key(ONEKG_REF_COL, ONEKG_ALT_COL).alias(
            PPP_INDEX_ALLELE_KEY_COL
        )
    )
    dbsnp_keyed = dbsnp.select(
        pl.col(GWASLAB_CHROM_COL).cast(pl.Int32),
        pl.col(GWASLAB_POS_COL).cast(pl.Int32),
        pl.col(GWASLAB_RSID_COL),
        unordered_allele_key(
            GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
        ).alias(PPP_INDEX_ALLELE_KEY_COL),
    ).unique(subset=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, PPP_INDEX_ALLELE_KEY_COL])

    return (
        filtered_keyed.join(
            dbsnp_keyed,
            on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, PPP_INDEX_ALLELE_KEY_COL],
            how="left",
        )
        .with_columns(
            pl.col(ONEKG_ALT_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
            pl.col(ONEKG_REF_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        )
        .select(MEMBERSHIP_COLUMNS)
    )


@frozen
class Common1kgMembershipTask(GeneratingTask):
    """
    Produce the normalized common-1kg membership list (CHR, POS, EA, NEA, rsID) by
    attaching nullable rsIDs to the filtered common 1kg variants.

    filtered_variants_task: a FilterCommon1kgVariantsTask output (CHR, POS, REF,
        ALT, AF).
    """

    meta: ReferenceFileMeta
    filtered_variants_task: Task

    @property
    def filtered_variants_meta(self) -> Meta:
        return self.filtered_variants_task.meta

    @property
    def deps(self) -> list[Task]:
        return [self.filtered_variants_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        filtered_asset = fetch(self.filtered_variants_task.asset_id)
        filtered = (
            scan_dataframe_asset(
                filtered_asset,
                meta=self.filtered_variants_meta,
                parquet_backend="polars",
            )
            .to_native()
            .collect()
        )
        dbsnp_path = gwaslab_download_ref_if_missing(
            GWASLAB_1KG_DBSNP151_HG38_AUTO_NAME
        )
        dbsnp = pl.read_csv(dbsnp_path, separator=_DBSNP_SEPARATOR)
        membership = attach_rsid(filtered, dbsnp)
        out_path = scratch_dir / "common_1kg_membership.parquet"
        write_byte_stream_split_parquet(membership, out_path, float_columns=[])
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, filtered_variants_task: Task, asset_id: str
    ) -> "Common1kgMembershipTask":
        return cls(
            meta=ReferenceFileMeta(
                group="ppp_variant_index_membership",
                sub_group="common_1kg",
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename="common_1kg_membership",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            filtered_variants_task=filtered_variants_task,
        )
