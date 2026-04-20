"""
Task to convert a JSON file to parquet using duckdb.
Implemented by claude
"""

from pathlib import Path

import duckdb
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


@frozen
class JsonToParquetTask(Task):
    """
    Task to convert a JSON file to a parquet file using duckdb.
    """

    meta: Meta
    json_task: Task
    target_compression: str = "zstd"

    @property
    def deps(self) -> list["Task"]:
        return [self.json_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.json_task.meta.asset_id)
        assert isinstance(source_asset, FileAsset)
        source_path = source_asset.path
        out_path = scratch_dir / "output_file.parquet"
        sql_command = f"""
            COPY (SELECT * FROM read_json_auto('{source_path}', maximum_object_size=1073741824))
            TO '{out_path}' (FORMAT 'PARQUET', CODEC '{self.target_compression}');
        """
        logger.info("running sql command", sql_command=sql_command)
        duckdb.sql(sql_command)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        json_task: Task,
        asset_id: str,
        target_compression: str = "zstd",
    ) -> "JsonToParquetTask":
        source_meta = json_task.meta
        if isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=source_meta.sub_folder,
                id=AssetId(asset_id),
                extension=f".parquet.{target_compression}",
                read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
            )
            return cls(
                meta=meta,
                json_task=json_task,
                target_compression=target_compression,
            )
        raise ValueError(f"Unsupported source meta type: {type(source_meta)}")
