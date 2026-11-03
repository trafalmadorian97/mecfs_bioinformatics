"""
Download a file, possibly verifying it using a hash.

This Task is the main way GWAS summary statistics are added to the build system.
"""

from pathlib import Path
from urllib.parse import urlsplit

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

_DROPBOX_PREVIEW_HOSTS = ("dropbox.com", "www.dropbox.com")


@frozen(slots=True)
class DownloadFileTask(Task):
    """
    connections optionally requests a number of simultaneous connections for the download.
    This can greatly speed up downloads from slow servers such as the GWAS Catalog.

    Dropbox links must use dl.dropboxusercontent.com. www.dropbox.com can answer dl=1 links
    with an HTML file preview page and HTTP status 200, depending on the client.
    """

    meta: Meta
    url: str
    md5_hash: str | None
    connections: int | None = None

    def __attrs_post_init__(self) -> None:
        assert self.connections is None or self.connections >= 1, (
            f"connections must be positive, got {self.connections}"
        )
        assert urlsplit(self.url).hostname not in _DROPBOX_PREVIEW_HOSTS, (
            f"Use dl.dropboxusercontent.com instead of www.dropbox.com: {self.url}"
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
