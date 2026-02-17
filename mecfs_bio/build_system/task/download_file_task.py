"""
Download a file, possibly verifying it using a hash.

This Task is the main GWAS summary statistics are added to the build system.
"""

from pathlib import Path

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


@frozen
class DownloadFileTask(Task):
    _meta: Meta
    _url: str
    _md5_hash: str | None

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        target = scratch_dir / self.meta.asset_id
        wf.download_from_url(url=self._url, local_path=target, md5_hash=self._md5_hash)
        return FileAsset(target)
