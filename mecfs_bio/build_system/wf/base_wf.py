from abc import ABC
from pathlib import Path

import py3_wget

from mecfs_bio.util.download.robust_download import robust_download_with_curl
from mecfs_bio.util.download.verify import verify_hash


class WF(ABC):
    """
    An interface to the external world.
    Currently only used for downloading files.
    May be extended in the future .
    """

    def download_from_url(
        self, url: str, local_path: Path, md5_hash: str | None
    ) -> None:
        pass


class SimpleWF(WF):
    def download_from_url(
        self, url: str, local_path: Path, md5_hash: str | None
    ) -> None:
        # py3_wget attempts to read the whole file to check md5.  doesn't work for large files.  So implement my own check
        py3_wget.download_file(
            url=url,
            output_path=local_path,
        )
        verify_hash(
            downloaded_file=local_path,
            expected_hash=md5_hash,
        )


class RobustDownloadWF(WF):
    def download_from_url(
        self, url: str, local_path: Path, md5_hash: str | None
    ) -> None:
        robust_download_with_curl(
            md5sum=md5_hash,
            url=url,
            dest=local_path,
        )
