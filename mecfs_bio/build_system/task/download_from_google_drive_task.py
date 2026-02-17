"""
Some GWAS summary statistics are stored on Google Drive.  Use this task to access them.
"""

from pathlib import Path

import gdown
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class DownloadFromGoogleDriveTask(Task):
    _meta: Meta
    file_id: str

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        out_path = scratch_dir / self.file_id
        gdown.download(id=self.file_id, output=str(out_path))
        return FileAsset(out_path)
