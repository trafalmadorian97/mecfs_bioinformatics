"""
Download a file, possibly verifying it using a hash.

This Task is the main way GWAS summary statistics are added to the build system.
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


@frozen(slots=True)
class DownloadFileTask(Task):
    """
    connections optionally requests a number of simultaneous connections for the download.
    This can greatly speed up downloads from slow servers such as the GWAS Catalog.
    """

    meta: Meta
    url: str
    md5_hash: str | None
    connections: int | None = None

    def __attrs_post_init__(self) -> None:
        assert self.connections is None or self.connections >= 1, (
            f"connections must be positive, got {self.connections}"
        )

    @property
    def deps(self) -> list["Task"]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        target = scratch_dir / self.meta.asset_id
        wf.download_from_url(
            url=self.url,
            local_path=target,
            md5_hash=self.md5_hash,
            requested_connections=self.connections,
        )
        return FileAsset(target)
