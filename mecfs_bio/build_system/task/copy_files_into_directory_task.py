import shutil
from pathlib import Path

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


@frozen
class CopySource:
    task: Task
    suffix: str | None = None
    name_in_dir: str | None = None


@frozen
class CopyFilesIntoDirectoryTask(Task):
    sources: list[CopySource]
    meta: Meta

    @property
    def deps(self) -> list["Task"]:
        return [item.task for item in self.sources]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        for source in self.sources:
            nm = source.name_in_dir or source.task.asset_id
            nm = nm if source.suffix is None else (nm + source.suffix)

            asset = fetch(source.task.asset_id)
            if isinstance(asset, FileAsset):
                logger.debug(
                    f"Copying file {source.task.asset_id} from {asset.path} to {(scratch_dir / nm,)}"
                )
                shutil.copy2(asset.path, scratch_dir / nm)
            elif isinstance(asset, DirectoryAsset):
                logger.debug(
                    f"Copying dir {source.task.asset_id} from {asset.path} to {(scratch_dir / nm,)}"
                )
                shutil.copytree(asset.path, scratch_dir / nm)
            else:
                raise ValueError(f"asset {asset} has unknown type")
        return DirectoryAsset(scratch_dir)
