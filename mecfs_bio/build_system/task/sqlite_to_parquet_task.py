"""
Task to query a SQLite database with duckdb and write the result as parquet.
Initial version written by Claude
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

_DB_ALIAS = "_src"


@frozen
class SqliteToParquetTask(Task):
    """
    Attaches a SQLite database as alias '_src', executes `query` against it,
    and writes the result to parquet.

    The `query` field must be a SELECT statement whose table references are
    prefixed with the alias '_src.' (e.g. 'SELECT * FROM _src.my_table').
    """

    meta: Meta
    source_task: Task
    query: str
    target_compression: str = "zstd"

    @property
    def deps(self) -> list["Task"]:
        return [self.source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_task.meta.asset_id)
        assert isinstance(source_asset, FileAsset)
        source_path = source_asset.path
        out_path = scratch_dir / "output_file.parquet"

        conn = duckdb.connect()
        conn.execute(f"ATTACH '{source_path}' AS {_DB_ALIAS} (TYPE SQLITE, READ_ONLY)")
        copy_sql = f"""
            COPY ({self.query})
            TO '{out_path}' (FORMAT 'PARQUET', CODEC '{self.target_compression}');
        """
        logger.info("running sql command", sql_command=copy_sql)
        conn.execute(copy_sql)
        conn.close()
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        source_task: Task,
        asset_id: str,
        query: str,
        target_compression: str = "zstd",
    ) -> "SqliteToParquetTask":
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
                query=query,
                target_compression=target_compression,
            )
        raise ValueError(f"Unsupported source meta type: {type(source_meta)}")
