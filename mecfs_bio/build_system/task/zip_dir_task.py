import shutil
from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_archive_meta import ResultArchiveMeta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class ZipDirTask(Task):
    source_dir_task: Task
    meta: Meta

    @property
    def deps(self) -> list["Task"]:
        return [self.source_dir_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_dir_task.asset_id)
        assert isinstance(source_asset, DirectoryAsset)
        source_path = source_asset.path
        shutil.make_archive(str(scratch_dir / "out"), "zip", source_path)
        return FileAsset(scratch_dir / "out.zip")

    @classmethod
    def create(
        cls,
        source_dir_task: Task,
        asset_id: str,
    ) -> "ZipDirTask":
        source_meta = source_dir_task.meta
        if isinstance(source_meta, ResultDirectoryMeta):
            meta = ResultArchiveMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                extension=".zip",
            )
        else:
            raise ValueError("Unknown source meta type")

        return ZipDirTask(
            source_dir_task=source_dir_task,
            meta=meta,
        )
