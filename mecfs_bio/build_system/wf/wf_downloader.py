"""
URL download capability
"""

from abc import ABC, abstractmethod
from pathlib import Path

import py3_wget

from mecfs_bio.util.download.robust_download import robust_download_with_aria
from mecfs_bio.util.download.verify import verify_hash


class WFDownloader(ABC):
    """Downloads the file at a URL to a local path, verifying its md5 hash."""

    @abstractmethod
    def download(self, url: str, local_path: Path, md5_hash: str | None) -> None:
        """Download the file at url into local_path and verify its md5 hash."""


class SimpleWFDownloader(WFDownloader):
    """Downloads via py3_wget, then verifies the hash separately."""

    def download(self, url: str, local_path: Path, md5_hash: str | None) -> None:
        # py3_wget attempts to read the whole file to check md5.  doesn't work for large files.  So implement my own check
        py3_wget.download_file(
            url=url,
            output_path=local_path,
        )
        verify_hash(
            downloaded_file=local_path,
            expected_hash=md5_hash,
        )


class RobustWFDownloader(WFDownloader):
    """Downloads via aria2 with retries, suitable for large files."""

    def download(self, url: str, local_path: Path, md5_hash: str | None) -> None:
        robust_download_with_aria(
            md5sum=md5_hash,
            url=url,
            dest=local_path,
        )
