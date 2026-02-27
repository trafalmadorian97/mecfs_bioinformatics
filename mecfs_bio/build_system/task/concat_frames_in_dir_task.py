"""
Concatenate dataframes in a DirectoryAsset to create a single FileAsset.
"""

from pathlib import Path, PurePath

import pandas as pd
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.getLogger()


@frozen
class ConcatFramesInDirTask(Task):
    """
    Task to concatenate dataframes in a directory matching a glob pattern

    Example use case: when GWAS summary statistics are split over multiple files.
    """

    _meta: Meta
    source_dir_task: Task
    path_glob: str
    read_spec_for_frames: DataFrameReadSpec

    @property
    def src_id(self) -> AssetId:
        return self.source_dir_task.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [self.source_dir_task]

    @property
    def meta(self) -> Meta:
        return self._meta

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        dir_asset = fetch(self.src_id)
        assert isinstance(dir_asset, DirectoryAsset)
        file_paths = sorted(dir_asset.path.rglob(self.path_glob))
        assert len(file_paths) > 0
        frames = []
        for pth in file_paths:
            frames.append(
                scan_dataframe(path=pth, spec=self.read_spec_for_frames)
                .collect()
                .to_pandas()
            )
            logger.debug(f"loaded frame from {pth}. shape: {frames[-1].shape}")
        result: pd.DataFrame = pd.concat(frames, axis=0)
        out_path = scratch_dir / "combined_dataframe.parquet"
        result.to_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_dir_task: Task,
        path_glob: str,
        read_spec_for_frames: DataFrameReadSpec,
    ):
        meta: Meta
        source_meta: Meta = source_dir_task.meta
        if isinstance(source_meta, ReferenceDataDirectoryMeta):
            meta = ReferenceFileMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=source_meta.sub_folder,
                id=AssetId(asset_id),
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            )
        elif isinstance(source_meta, ProcessedGwasDataDirectoryMeta):
            meta = FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=PurePath("processed"),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            )
        else:
            raise NotImplementedError(
                f"Handlers for meta {source_meta} are not implemented"
            )
        return cls(
            meta=meta,
            source_dir_task=source_dir_task,
            path_glob=path_glob,
            read_spec_for_frames=read_spec_for_frames,
        )
