from pathlib import Path
from typing import Callable, Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.download_file_task import calc_md5_checksum
from mecfs_bio.build_system.wf.base_wf import WF


import structlog

from mecfs_bio.util.external_util.bcftools import build_bcftools_index_command
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()


@frozen
class DownloadEntry:
    """
    A file to download into the directory,
    possibly including a command line operation to call on the downloaded file
    """
    url: str
    filename:str
    md5hash: str
    post_download_command: Callable[[Path],list[str]] | None
    @classmethod
    def create_with_bcftools_index(cls,
                                   url:str,
                                   filename:str,
                                   md5hash:str):
        return cls(
            url=url,
            filename=filename,
            md5hash=md5hash,
            post_download_command=build_bcftools_index_command
        )


@frozen
class DownloadFilesIntoDirectoryTask(Task):
    """
    A directory consisting of one or more downloaded files
    """

    _meta: Meta
    entries: Sequence[DownloadEntry]

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        for entry in self.entries:
            target = scratch_dir / entry.filename
            logger.debug(f"Downloading {entry}")
            wf.download_from_url(entry.url, target)
            logger.debug("Verifying MD5 hash of downloaded file...")
            hash_of_downloaded_file = calc_md5_checksum(target)
            assert hash_of_downloaded_file == entry.md5hash, (
                f"Expected Hash {hash_of_downloaded_file} to be equal to {entry.md5hash}"
            )
            logger.debug("Hash verified.")
            if entry.post_download_command is not None:
                post_download_command = entry.post_download_command(target)
                logger.debug(f"Running post_download_command: {" ".join(post_download_command)} \n")
                execute_command(post_download_command)
        return DirectoryAsset(scratch_dir)



