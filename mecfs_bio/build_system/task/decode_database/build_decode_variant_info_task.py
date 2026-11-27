"""
deCODE's variant_info_iceland table (allele frequencies per allele): the raw download, a sorted
parquet built from it, and the clean subset of that parquet.

DownloadDecodeVariantInfoTask keeps the source as published (gzipped tab-separated text,
Chrom Pos A1 ImpFreqA1 ImpInfo), in its original row order. The object carries no md5
metadata; its ETag (a single-part upload) is its md5 and is checked after the download.

BuildDecodeVariantInfoTask: the source writes most rows several times; the parquet keeps each
distinct (CHR, POS, A1, ImpFreqA1, ImpInfo) row once, sorted by (CHR, POS) in small row groups
so that readers of one window can skip most of the file.

BuildDecodeCleanVariantInfoTask keeps the parquet's entries at positions where it has exactly
one entry. A deCODE row is clean when its (CHR, POS, effect allele) is one of these entries
(keep_clean_rows): variant_info then lists nothing at the row's position but the row's own
effect allele. Only clean rows are two-allele tests whose effect allele is contrasted with the
reference allele alone. variant_info lists a two-allele record's ALT only, but a multi-allelic
record's every allele (REF and "*" included), and the entries of records overlapping one
position are interleaved; at such positions deCODE's rows are one allele against the pooled
rest, which a trait GWAS's REF/ALT row does not measure. A repeated allele with different
frequencies (two records) is not clean either. On the AP4M1 file the rule keeps 93.2% of rows
with MAF >= 1% that are SNVs, 99.5% of SNVs alone at their position, and 23.9% of indels with
MAF >= 1%.
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
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

_KEY_CHR_COL = "__chr__"
_KEY_POS_COL = "__pos__"
_KEY_ALLELE_COL = "__allele__"
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


@frozen(slots=True)
class BuildDecodeCleanVariantInfoTask(GeneratingTask):
    """variant_info_task: the BuildDecodeVariantInfoTask of the deduplicated parquet."""

    meta: ReferenceFileMeta
    variant_info_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.variant_info_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.variant_info_task.asset_id)
        assert isinstance(source, FileAsset)
        out_path = scratch_dir / "variant_info_iceland_clean.parquet"
        write_clean_variant_info_parquet(source.path, out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, asset_id: str, variant_info_task: Task
    ) -> "BuildDecodeCleanVariantInfoTask":
        return cls(
            meta=ReferenceFileMeta(
                group="decode_somascan",
                sub_group="variant_info",
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename="variant_info_iceland_clean",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            variant_info_task=variant_info_task,
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


def write_clean_variant_info_parquet(source: Path, out_path: Path) -> None:
    """The entries of the deduplicated variant_info parquet at source at single-entry positions."""
    (
        pl.read_parquet(source)
        .filter(pl.len().over(VI_CHR_COL, VI_POS_COL) == 1)
        .write_parquet(
            out_path, compression="zstd", row_group_size=VARIANT_INFO_ROW_GROUP_SIZE
        )
    )


def keep_clean_rows(
    rows: pl.DataFrame, clean_variant_info: pl.DataFrame
) -> pl.DataFrame:
    """
    rows: CHR (chromosome code), POS and EA, among others.
    clean_variant_info: entries (CHR, POS, A1, ...) of the clean table at least at the rows'
        positions.
    Returns the clean rows, in their order.
    """
    keys = clean_variant_info.select(
        pl.col(VI_CHR_COL).cast(pl.Int8).alias(_KEY_CHR_COL),
        pl.col(VI_POS_COL).cast(pl.Int64).alias(_KEY_POS_COL),
        pl.col(VI_A1_COL).alias(_KEY_ALLELE_COL),
    )
    return (
        rows.with_columns(
            pl.col(GWASLAB_CHROM_COL).cast(pl.Int8).alias(_KEY_CHR_COL),
            pl.col(GWASLAB_POS_COL).cast(pl.Int64).alias(_KEY_POS_COL),
            pl.col(GWASLAB_EFFECT_ALLELE_COL).alias(_KEY_ALLELE_COL),
        )
        .join(
            keys,
            on=[_KEY_CHR_COL, _KEY_POS_COL, _KEY_ALLELE_COL],
            how="semi",
            maintain_order="left",
        )
        .drop(_KEY_CHR_COL, _KEY_POS_COL, _KEY_ALLELE_COL)
    )
