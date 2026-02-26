"""
Task to combine the results of multiple Tasks, each of which produces a dataframe.
"""
from pathlib import Path
from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat, OutFormat, CSVOutFormat
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class ConcatFramesTask(Task):
    _meta: Meta
    frames_tasks: Sequence[Task]
    out_format:  OutFormat

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return list(self.frames_tasks)

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        frames = []
        for task in self.frames_tasks:
            asset  = fetch(task.asset_id)
            frames.append(scan_dataframe_asset(asset, meta=task.meta))
        result = narwhals.concat(frames,how="vertical")
        out_path = scratch_dir / f"{self.meta.asset_id}"
        if isinstance(self.out_format, CSVOutFormat):
            result.collect().to_pandas().to_csv(
                out_path, index=False, sep=self.out_format.sep
            )
        elif isinstance(self.out_format, ParquetOutFormat):
            result.sink_parquet(out_path)
        return FileAsset(out_path)
