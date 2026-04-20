"""
Task to unpack a parquet file with a MAP(VARCHAR, STRUCT) column into a flat table.
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
class UnpackMapParquetTask(Task):
    """
    Reads a parquet whose column is MAP(VARCHAR, STRUCT(...)), unnests the map entries
    into one row per key, and expands the struct fields as individual columns.

    The map keys are emitted as `name_column`; all struct fields follow as columns.
    """

    meta: Meta
    source_task: Task
    map_column: str = "json"
    name_column: str = "name"
    target_compression: str = "zstd"
    def __attrs_post_init__(self):
        assert isinstance(self.source_task.meta.read_spec, DataFrameReadSpec)
        assert isinstance(self.source_task.meta.read_spec.format, DataFrameParquetFormat)

    @property
    def deps(self) -> list["Task"]:
        return [self.source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_task.meta.asset_id)
        assert isinstance(source_asset, FileAsset)
        source_path = source_asset.path
        out_path = scratch_dir / "output_file.parquet"
        sql_command = f"""
            COPY (
                SELECT {self.name_column}, val.*
                FROM (
                    SELECT
                        unnest(map_keys({self.map_column})) AS {self.name_column},
                        unnest(map_values({self.map_column})) AS val
                    FROM read_parquet('{source_path}')
                )
            )
            TO '{out_path}' (FORMAT 'PARQUET', CODEC '{self.target_compression}');
        """
        logger.info("running sql command", sql_command=sql_command)
        duckdb.sql(sql_command)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        source_task: Task,
        asset_id: str,
        map_column: str = "json",
        name_column: str = "name",
        target_compression: str = "zstd",
    ) -> "UnpackMapParquetTask":
        source_meta = source_task.meta
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
                source_task=source_task,
                map_column=map_column,
                name_column=name_column,
                target_compression=target_compression,
            )
        raise ValueError(f"Unsupported source meta type: {type(source_meta)}")
