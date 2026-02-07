import hashlib
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
        wf.download_from_url(url=self._url, local_path=target)
        verify_hash(target, expected_hash=self._md5_hash)
        return FileAsset(target)


def verify_hash(downloaded_file: Path, expected_hash: str | None):
    if expected_hash is None:
        return
    logger.debug("Verifying MD5 hash of downloaded file...")
    hash_of_downloaded_file = calc_md5_checksum(downloaded_file)
    if hash_of_downloaded_file == expected_hash:
        logger.debug("Hash verified.")
        return
    head_file(downloaded_file)
    raise AssertionError(
        f"Expected has {hash_of_downloaded_file} of file {downloaded_file} to be equal to {expected_hash}"
    )


def head_file(filename: Path, n=10):
    """Prints the first n lines of a file, like the Unix head command."""
    logger.debug(f"first {n} lines of file {filename}:")
    try:
        with open(filename) as f:
            for i, line in enumerate(f):
                if i >= n:
                    break
                logger.debug(line.rstrip("\n"))  # rstrip to avoid double newlines
    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")


def calc_md5_checksum(filepath: Path, chunk_size: int = 8192) -> str:
    hasher = hashlib.md5()
    with open(filepath, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break  # End of file
            hasher.update(chunk)
    return hasher.hexdigest()
