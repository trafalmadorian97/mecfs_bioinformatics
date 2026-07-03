"""
Task that converts an arbitrary-whitespace-separated (optionally gzipped) text
table into parquet, streaming in bounded-memory chunks.

Polars cannot read arbitrary-whitespace-separated files, so the read_spec machinery
falls back to a single eager pandas read for DataFrameWhiteSpaceSepTextFormat (see
mecfs_bio/build_system/meta/read_spec/read_dataframe.py). For a large meta-analysis
such as the deCODE summary statistics (tens of millions of rows, object-dtype string
columns), that materialises the whole table in pandas and then copies it into polars,
which exhausts memory.

Converting the file to parquet once, up front, lets every downstream task scan it
lazily with projection pushdown. The conversion itself reads the source in row
chunks and appends each chunk to the parquet file, so peak memory is bounded by
chunk_size rather than the size of the whole table. The chunked read uses the same
pandas whitespace parsing (sep=r"\\s+") as the read_spec machinery, so the parsed
result is identical to reading the file directly.
"""

from pathlib import Path, PurePath
from typing import Mapping

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


@frozen
class WhitespaceSepTextToParquetTask(Task):
    meta: Meta
    source_task: Task
    chunk_size: int = 2_000_000
    compression: str = "zstd"
    dtype: Mapping[str, str] | None = None

    def __attrs_post_init__(self) -> None:
        # Shift-left: reject a source that is not arbitrary-whitespace-separated text
        # at construction time rather than waiting until execute.
        self._whitespace_format()

    @property
    def deps(self) -> list["Task"]:
        return [self.source_task]

    @property
    def _source_id(self) -> AssetId:
        return self.source_task.asset_id

    @property
    def _source_meta(self) -> FileMeta:
        meta = self.source_task.meta
        assert isinstance(meta, FileMeta)
        return meta

    def _whitespace_format(self) -> DataFrameWhiteSpaceSepTextFormat:
        read_spec = self._source_meta.read_spec
        assert read_spec is not None, "source asset has no read_spec"
        fmt = read_spec.format
        assert isinstance(fmt, DataFrameWhiteSpaceSepTextFormat), (
            "source asset must be arbitrary-whitespace-separated text"
        )
        return fmt

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self._source_id)
        assert isinstance(source_asset, FileAsset)
        fmt = self._whitespace_format()

        extra_options: dict = {}
        if fmt.col_names is not None:
            extra_options["names"] = fmt.col_names
        if self.dtype is not None:
            extra_options["dtype"] = self.dtype

        out_path = scratch_dir / "output_file.parquet"
        writer: pq.ParquetWriter | None = None
        schema: pa.Schema | None = None
        rows = 0
        # sep=r"\s+" matches read_dataframe.scan_dataframe; pandas infers gzip from the
        # path extension. The first chunk fixes the schema; later chunks are cast to it
        # (safe=False) so an all-null column that pandas widens later cannot drift.
        reader = pd.read_csv(
            source_asset.path,
            sep=r"\s+",
            comment=fmt.comment_code,
            chunksize=self.chunk_size,
            **extra_options,
        )
        for chunk in reader:
            if schema is None:
                table = pa.Table.from_pandas(chunk, preserve_index=False)
                schema = table.schema
                writer = pq.ParquetWriter(
                    out_path, schema, compression=self.compression
                )
            else:
                table = pa.Table.from_pandas(
                    chunk, schema=schema, preserve_index=False, safe=False
                )
            assert writer is not None
            writer.write_table(table)
            rows += chunk.shape[0]
        assert writer is not None, "source table is empty"
        writer.close()
        logger.info(f"Converted {rows} rows to parquet at {out_path}")
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        source_task: Task,
        asset_id: str,
        sub_dir: str | PurePath = PurePath("processed"),
        chunk_size: int = 2_000_000,
        compression: str = "zstd",
        dtype: Mapping[str, str] | None = None,
    ) -> "WhitespaceSepTextToParquetTask":
        source_meta = source_task.meta
        assert isinstance(
            source_meta, (GWASSummaryDataFileMeta, FilteredGWASDataMeta)
        ), "source task must carry GWAS summary metadata"
        meta = FilteredGWASDataMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=sub_dir,
            read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
        )
        return cls(
            meta=meta,
            source_task=source_task,
            chunk_size=chunk_size,
            compression=compression,
            dtype=dtype,
        )
