import shutil
from pathlib import Path, PurePath

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.procesed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class CopyFileFromDirectoryTask(Task):
    _meta: Meta
    source_directory_task: Task
    path_inside_directory: PurePath

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.source_directory_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        dir_asset = fetch(self.source_directory_task.asset_id)
        assert isinstance(dir_asset, DirectoryAsset)
        src_path = dir_asset.path / self.path_inside_directory
        out_path = scratch_dir / self.meta.asset_id
        shutil.copyfile(src_path, out_path)
        return FileAsset(out_path)

    @classmethod
    def create_result_table(
        cls,
        asset_id: str,
        source_directory_task: Task,
        path_inside_directory: PurePath,
        extension: str,
        read_spec: DataFrameReadSpec | None,
    ):
        source_meta = source_directory_task.meta
        if isinstance(source_meta, ProcessedGwasDataDirectoryMeta):
            meta = ResultTableMeta(
                asset_id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                extension=extension,
                read_spec=read_spec,
            )
        else:
            raise ValueError("Unknown source meta")
        return cls(
            meta=meta,
            source_directory_task=source_directory_task,
            path_inside_directory=path_inside_directory,
        )
