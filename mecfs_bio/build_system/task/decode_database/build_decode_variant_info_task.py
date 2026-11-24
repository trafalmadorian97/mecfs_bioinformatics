"""
deCODE's variant_info_iceland table (allele frequencies per allele): the raw download, and a
sorted parquet built from it.

DownloadDecodeVariantInfoTask keeps the source as published (gzipped tab-separated text,
Chrom Pos A1 ImpFreqA1 ImpInfo), in its original row order. The object carries no md5
metadata; its ETag (a single-part upload) is its md5 and is checked after the download.

BuildDecodeVariantInfoTask: the source writes most rows several times; the parquet keeps each
distinct (CHR, POS, A1, ImpFreqA1, ImpInfo) row once, sorted by (CHR, POS) in small row groups
so that readers of one window can skip most of the file.
"""

from pathlib import Path, PurePath

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code_expr,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.build_system.wf.s3_access import S3ObjectRef
from mecfs_bio.constants.decode_somascan_constants import (
    VI_A1_COL,
    VI_CHR_COL,
    VI_FREQ_COL,
    VI_INFO_COL,
    VI_POS_COL,
)

# Column names in the source file.
_SOURCE_CHROM_COL = "Chrom"
_SOURCE_POS_COL = "Pos"
VARIANT_INFO_ROW_GROUP_SIZE = 65_536
_PRESIGN_EXPIRES_S = 6 * 3600


@frozen(slots=True)
class DownloadDecodeVariantInfoTask(GeneratingTask):
    """
    ref: the variant_info_iceland.txt.gz object.
    n_connections: parallel range requests for the download.
    """

    meta: ReferenceFileMeta
    ref: S3ObjectRef
    n_connections: int = 128

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        download = scratch_dir / "variant_info_iceland.txt.gz"
        _download_etag_verified(self.ref, download, self.n_connections, wf)
        return FileAsset(download)

    @classmethod
    def create(cls, asset_id: str, ref: S3ObjectRef) -> "DownloadDecodeVariantInfoTask":
        return cls(
            meta=ReferenceFileMeta(
                group="decode_somascan",
                sub_group="variant_info",
                sub_folder=PurePath("raw"),
                id=AssetId(asset_id),
                filename="variant_info_iceland",
                extension=".txt.gz",
                read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
            ),
            ref=ref,
        )


@frozen(slots=True)
class BuildDecodeVariantInfoTask(GeneratingTask):
    """source_task: the DownloadDecodeVariantInfoTask of the raw text."""

    meta: ReferenceFileMeta
    source_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.source_task.asset_id)
        assert isinstance(source, FileAsset)
        out_path = scratch_dir / "variant_info_iceland.parquet"
        write_variant_info_parquet(source.path, out_path)
        return FileAsset(out_path)

    @classmethod
    def create(cls, asset_id: str, source_task: Task) -> "BuildDecodeVariantInfoTask":
        return cls(
            meta=ReferenceFileMeta(
                group="decode_somascan",
                sub_group="variant_info",
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename="variant_info_iceland",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            source_task=source_task,
        )


def _download_etag_verified(
    ref: S3ObjectRef, local_path: Path, n_connections: int, wf: WF
) -> None:
    info = wf.s3_access.head(ref)
    assert info.md5 is None, (
        f"{ref.key} now carries md5 metadata; check it instead of the ETag"
    )
    wf.download_url_parallel(
        url=wf.s3_access.presigned_url(ref, _PRESIGN_EXPIRES_S),
        size=info.size,
        local_path=local_path,
        md5_hash=info.etag,
        n_connections=n_connections,
    )


def write_variant_info_parquet(source: Path, out_path: Path) -> None:
    """Deduplicate, code, cast and sort the variant_info text at source into out_path."""
    (
        pl.scan_csv(source, separator="\t")
        .select(
            chrom_code_expr(_SOURCE_CHROM_COL).alias(VI_CHR_COL),
            pl.col(_SOURCE_POS_COL).cast(pl.Int64).alias(VI_POS_COL),
            pl.col(VI_A1_COL).cast(pl.String),
            pl.col(VI_FREQ_COL).cast(pl.Float64),
            pl.col(VI_INFO_COL).cast(pl.Float32),
        )
        .unique()
        .sort(VI_CHR_COL, VI_POS_COL)
        .sink_parquet(
            out_path, compression="zstd", row_group_size=VARIANT_INFO_ROW_GROUP_SIZE
        )
    )
