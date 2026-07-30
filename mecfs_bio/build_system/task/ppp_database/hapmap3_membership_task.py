"""
Normalized HapMap3 membership list for the PPP variant index (HapMap3 mode).

Reads the HapMap3 hg38 snplist bundled inside the gwaslab package -- guaranteed
present because gwaslab is a pinned pyproject.toml dependency, so no download or
rehosting is needed -- and normalizes it to the membership contract
(CHR, POS[hg38], EA, NEA, rsID). Allele orientation here is arbitrary (the index
adopts the template protein's orientation); only the allele set is used for matching.
"""

import os
from pathlib import Path, PurePath

import gwaslab as gl
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.ppp_database.byte_stream_split_parquet import (
    write_byte_stream_split_parquet,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_HAPMAP3_HG38_SNPLIST_RELPATH,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)

# Raw column names in the gwaslab HapMap3 snplist.
_HM3_RSID_COL = "rsid"
_HM3_A1_COL = "A1"
_HM3_A2_COL = "A2"
_HM3_CHROM_COL = "#CHROM"
_HM3_POS_COL = "POS"

MEMBERSHIP_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
]


def bundled_hapmap3_hg38_snplist_path() -> Path:
    """Path to the HapMap3 hg38 snplist bundled in the installed gwaslab package."""
    return Path(os.path.dirname(gl.__file__)) / GWASLAB_HAPMAP3_HG38_SNPLIST_RELPATH


def normalize_hapmap3_snplist(snplist_path: Path) -> pl.DataFrame:
    """Read the gwaslab HapMap3 snplist and normalize to the membership contract."""
    return (
        pl.read_csv(snplist_path, separator="\t")
        .select(
            pl.col(_HM3_CHROM_COL).cast(pl.Int32).alias(GWASLAB_CHROM_COL),
            pl.col(_HM3_POS_COL).cast(pl.Int32).alias(GWASLAB_POS_COL),
            pl.col(_HM3_A1_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
            pl.col(_HM3_A2_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
            pl.col(_HM3_RSID_COL).alias(GWASLAB_RSID_COL),
        )
        .select(MEMBERSHIP_COLUMNS)
    )


@frozen
class Hapmap3MembershipTask(GeneratingTask):
    """
    Produce the normalized HapMap3 membership list (CHR, POS, EA, NEA, rsID) from
    the gwaslab-bundled hg38 snplist.
    """

    meta: ReferenceFileMeta

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        membership = normalize_hapmap3_snplist(bundled_hapmap3_hg38_snplist_path())
        out_path = scratch_dir / "hapmap3_membership.parquet"
        write_byte_stream_split_parquet(membership, out_path, float_columns=[])
        return FileAsset(out_path)

    @classmethod
    def create(cls, asset_id: str) -> "Hapmap3MembershipTask":
        return cls(
            meta=ReferenceFileMeta(
                group="ppp_variant_index_membership",
                sub_group="hapmap3",
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename="hapmap3_membership",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            )
        )
