"""
Task to transform a dataframe asset using a DataProcessingPipes.
"""

from pathlib import Path, PurePath
from typing import Literal, Sequence

import pyarrow
import pyarrow.parquet
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
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    ValidBackend,
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.wf.base_wf import WF

ParquetCompression = Literal["snappy", "zstd", "gzip", "brotli", "lz4", "none"]


@frozen
class ParquetWriteOptions:
    """How to encode a parquet output, when the defaults are not good enough.

    Supplying these switches the write from polars sink_parquet to the pyarrow
    writer, which exposes per-column encoding control that sink_parquet does
    not. The trade-off is that the pyarrow path must materialize the frame in
    memory, so prefer the default (no options) for large outputs.

    byte_stream_split_floats is the reason this exists. It transposes the bytes
    of each floating-point value so that like-significance bytes sit together,
    which a general-purpose compressor exploits far better than interleaved
    IEEE754 bytes. On a table of GWAS statistics this cuts the compressed size
    by roughly a third over plain encoding, which matters for tables shipped to
    the documentation site and downloaded by readers' browsers.
    """

    compression: ParquetCompression = "zstd"
    compression_level: int | None = None
    byte_stream_split_floats: bool = False


@frozen
class ParquetOutFormat:
    write_options: ParquetWriteOptions | None = None


@frozen
class CSVOutFormat:
    sep: str


OutFormat = ParquetOutFormat | CSVOutFormat


@frozen
class PipeDataFrameTask(Task):
    source_data_task: Task
    pipes: Sequence[DataProcessingPipe]
    meta: Meta
    out_format: OutFormat
    backend: ValidBackend = "ibis"

    def __attrs_post_init__(self):
        source_meta = self._source_meta
        assert isinstance(source_meta, FileMeta)
        read_spec = source_meta.read_spec
        assert read_spec is not None
        if isinstance(read_spec.format, DataFrameTextFormat):
            assert self.backend in ("polars",), "Can only read text data with polars"

    @property
    def deps(self) -> list["Task"]:
        return [self.source_data_task]

    @property
    def _source_id(self) -> AssetId:
        return self.source_data_task.asset_id

    @property
    def _source_meta(self) -> Meta:
        return self.source_data_task.meta

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        asset = fetch(self._source_id)
        df = scan_dataframe_asset(
            asset=asset, meta=self._source_meta, parquet_backend=self.backend
        )
        out_path = scratch_dir / "out_dataframe"
        for pipe in self.pipes:
            df = pipe.process(df)
        if isinstance(self.out_format, CSVOutFormat):
            df.collect().to_pandas().to_csv(
                out_path, index=False, sep=self.out_format.sep
            )
        elif isinstance(self.out_format, ParquetOutFormat):
            if self.out_format.write_options is None:
                df.sink_parquet(out_path)
            else:
                write_parquet_with_options(
                    table=df.collect().to_arrow(),
                    out_path=out_path,
                    options=self.out_format.write_options,
                )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        source_task: Task,
        asset_id: str,
        out_format: OutFormat,
        pipes: Sequence[DataProcessingPipe],
        backend: ValidBackend = "ibis",
    ) -> "PipeDataFrameTask":
        source_meta = source_task.meta
        extension, read_spec = get_extension_and_read_spec_from_format(
            out_format=out_format
        )
        meta: Meta
        if isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                extension=extension,
                read_spec=read_spec,
            )
        elif isinstance(source_meta, GWASSummaryDataFileMeta):
            meta = FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=PurePath("processed"),
                read_spec=read_spec,
            )
        elif isinstance(source_meta, ResultTableMeta):
            meta = ResultTableMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                extension=extension,
                read_spec=read_spec,
            )
        elif isinstance(source_meta, FilteredGWASDataMeta):
            meta = FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=source_meta.sub_dir,
                read_spec=read_spec,
            )
        else:
            raise ValueError("unknown source meta")

        return cls(
            source_data_task=source_task,
            pipes=list(pipes),
            meta=meta,
            out_format=out_format,
            backend=backend,
        )


def write_parquet_with_options(
    table: pyarrow.Table, out_path: Path, options: ParquetWriteOptions
) -> None:
    """Write a parquet file with explicit encoding control.

    Dictionary encoding takes precedence over BYTE_STREAM_SPLIT in the parquet
    writer: a column left dictionary-enabled is written as RLE_DICTIONARY and
    the requested split is silently dropped, producing a file byte-identical to
    one written without it. Dictionary encoding is therefore disabled on exactly
    the float columns and left on for the rest, where it is what makes
    low-cardinality string columns small.
    """
    float_cols = [
        field.name for field in table.schema if pyarrow.types.is_floating(field.type)
    ]
    other_cols = [
        field.name
        for field in table.schema
        if not pyarrow.types.is_floating(field.type)
    ]
    split = options.byte_stream_split_floats
    pyarrow.parquet.write_table(
        table,
        out_path,
        compression=options.compression,
        compression_level=options.compression_level,
        use_byte_stream_split=float_cols if split else False,
        use_dictionary=other_cols if split else True,
    )


def get_extension_and_read_spec_from_format(
    out_format: OutFormat,
) -> tuple[str, DataFrameReadSpec]:
    if isinstance(out_format, CSVOutFormat):
        read_spec = DataFrameReadSpec(DataFrameTextFormat(separator=out_format.sep))
        if out_format.sep == "\t":
            extension = ".tsv"
        elif out_format.sep == ",":
            extension = ".csv"
        else:
            raise ValueError("Unknown sep")
    elif isinstance(out_format, ParquetOutFormat):
        read_spec = DataFrameReadSpec(DataFrameParquetFormat())
        extension = ".parquet"
    else:
        raise ValueError(f"Unknown format {out_format}")
    return extension, read_spec
