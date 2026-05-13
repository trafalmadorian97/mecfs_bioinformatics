"""
Task to combine the results of multiple Tasks, each of which produces a dataframe.
"""

from pathlib import Path
from typing import Mapping, Sequence

import narwhals
from attrs import frozen
from frozendict import frozendict

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    CSVOutFormat,
    OutFormat,
    ParquetOutFormat,
    get_extension_and_read_spec_from_format,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class ConcatFramesTask(Task):
    """
    Task to concatenate multiple DataFrames, each produces by a separate task.
    """

    meta: Meta
    frames_tasks: Sequence[Task]
    out_format: OutFormat
    frames_pipes: Sequence[DataProcessingPipe] | None = None
    column_type_override: Mapping[str, narwhals.dtypes.DType] = frozendict()

    def __attrs_post_init__(self):
        if self.frames_pipes is not None:
            assert len(self.frames_pipes) == len(self.frames_tasks)

    @property
    def deps(self) -> list["Task"]:
        return list(self.frames_tasks)

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        frames = []
        for i, task in enumerate(self.frames_tasks):
            asset = fetch(task.asset_id)
            frame = scan_dataframe_asset(asset, meta=task.meta)
            if len(self.column_type_override) > 0:
                frame = frame.with_columns(
                    *[
                        narwhals.col(col).cast(t)
                        for col, t in self.column_type_override.items()
                    ]
                )
            if self.frames_pipes is not None:
                frame = self.frames_pipes[i].process(frame)
            frames.append(frame)
        result = narwhals.concat(frames, how="vertical")
        out_path = scratch_dir / f"{self.meta.asset_id}"
        if isinstance(self.out_format, CSVOutFormat):
            result.collect().to_pandas().to_csv(
                out_path, index=False, sep=self.out_format.sep
            )
        elif isinstance(self.out_format, ParquetOutFormat):
            result.sink_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        frames_tasks: Sequence[Task],
        out_format: OutFormat,
        override_trait: str | None = None,
        override_project: str | None = None,
        column_type_override: Mapping[str, narwhals.dtypes.DType] = frozendict(),
        frames_pipes: Sequence[DataProcessingPipe] | None = None,
    ):
        extension, spec = get_extension_and_read_spec_from_format(out_format)
        assert len(frames_tasks) > 0
        source_meta = frames_tasks[0].meta
        if isinstance(source_meta, ResultTableMeta):
            trait = override_trait if override_trait is not None else source_meta.trait
            project = (
                override_project
                if override_project is not None
                else source_meta.project
            )
            meta = ResultTableMeta(
                id=AssetId(asset_id),
                trait=trait,
                project=project,
                extension=extension,
                read_spec=spec,
                sub_dir=source_meta.sub_dir,
            )
        else:
            raise ValueError(f"Unknown meta: {source_meta}")

        return cls(
            meta=meta,
            frames_tasks=frames_tasks,
            out_format=out_format,
            column_type_override=column_type_override,
            frames_pipes=frames_pipes,
        )
