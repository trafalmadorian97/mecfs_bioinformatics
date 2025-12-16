import structlog

from mecfs_bio.build_system.meta.base_meta import FileMeta

logger = structlog.get_logger()
from pathlib import Path

import duckdb
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class CompressedCSVToParquetTask(Task):
    """
    Task for converting a gzipped CSV-type file to a parquet file
    Main use is for processing the SNP151 SNP database files
    """

    _meta: Meta
    csv_task: Task
    source_compression: str = "gzip"
    target_compression: str = "zstd"
    select_list: list[str] | None = None

    @property
    def _source_meta(self) -> FileMeta:
        meta = self.csv_task.meta
        assert isinstance(meta, FileMeta)
        return meta

    @property
    def _source_id(self) -> AssetId:
        return self.csv_task.meta.asset_id

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.csv_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self._source_id)
        assert isinstance(source_asset, FileAsset)
        source_path = source_asset.path
        read_spec = self._source_meta.read_spec()
        format = _get_format(read_spec)
        assert format.null_values is None
        delim = _get_sep(format)
        col_names = _get_column_names(format)
        out_path = scratch_dir / "output_file.parquet"
        name_list_str = "["
        for nm in col_names:
            name_list_str += f"'{nm}',"
        name_list_str += "]"
        select_list = "*" if self.select_list is None else ",".join(self.select_list)
        sql_command = f"""
                    COPY (SELECT {select_list}
        	        FROM read_csv('{source_path}',
        	        AUTO_DETECT=TRUE, 
        	        HEADER={format.has_header}, 
        	        NAMES = {name_list_str}, 
        	        delim = '{delim}', 
        	        """
        if format.comment_char is not None:
            sql_command += f"""comment = '{format.comment_char}',
            """

        sql_command += f"""compression={self.source_compression} ))
                    TO '{out_path}' (FORMAT 'PARQUET', CODEC '{self.target_compression}');
                    """

        logger.info(f"running sql command:\n {sql_command}")
        duckdb.sql(sql_command)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        csv_task: Task,
        asset_id: str,
        target_compression: str = "zstd",
        source_compression: str = "gzip",
        select_list: list[str] | None = None,
    ) -> "CompressedCSVToParquetTask":
        source_meta = csv_task.meta
        if isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=source_meta.sub_folder,
                asset_id=AssetId(asset_id),
                extension=f".parquet.{target_compression}",
                read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
            )
            return cls(
                meta=meta,
                csv_task=csv_task,
                target_compression=target_compression,
                source_compression=source_compression,
                select_list=select_list,
            )
        raise ValueError("Unknown Meta")


def _get_sep(format: DataFrameTextFormat) -> str:
    return format.separator


def _get_format(spec: DataFrameReadSpec | None) -> DataFrameTextFormat:
    assert spec is not None
    format = spec.format
    assert isinstance(format, DataFrameTextFormat)
    return format


def _get_column_names(format: DataFrameTextFormat) -> list[str]:
    cols = format.column_names
    assert cols is not None
    return cols
