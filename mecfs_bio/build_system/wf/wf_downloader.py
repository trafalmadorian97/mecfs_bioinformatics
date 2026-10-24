"""
URL download capability
"""

from abc import ABC, abstractmethod
from pathlib import Path

import py3_wget
import structlog

from mecfs_bio.util.download.robust_download import (
    ARIA2_MAX_CONNECTIONS_PER_SERVER,
    robust_download_with_aria,
)
from mecfs_bio.util.download.verify import verify_hash

logger = structlog.get_logger()


class WFDownloader(ABC):
    """Downloads the file at a URL to a local path, verifying its md5 hash."""

    @abstractmethod
    def download(
        self,
        url: str,
        local_path: Path,
        md5_hash: str | None,
        requested_connections: int | None,
    ) -> None:
        """
        Download the file at url into local_path and verify its md5 hash.

        requested_connections is a hint for the number of simultaneous connections to use.
        Implementations may ignore it.
        """


class SimpleWFDownloader(WFDownloader):
    """Downloads via py3_wget, then verifies the hash separately.  Ignores requested_connections."""

    def download(
        self,
        url: str,
        local_path: Path,
        md5_hash: str | None,
        requested_connections: int | None,
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


class RobustWFDownloader(WFDownloader):
    """
    Downloads via aria2 with retries, suitable for large files.

    Honors requested_connections, capped at ARIA2_MAX_CONNECTIONS_PER_SERVER.
    """

    def download(
        self,
        url: str,
        local_path: Path,
        md5_hash: str | None,
        requested_connections: int | None,
    ) -> None:
        robust_download_with_aria(
            md5sum=md5_hash,
            url=url,
            dest=local_path,
            num_simil=_aria_connections(requested_connections),
        )


def _aria_connections(requested_connections: int | None) -> int:
    if requested_connections is None:
        return 1
    if requested_connections > ARIA2_MAX_CONNECTIONS_PER_SERVER:
        logger.warning(
            f"Requested {requested_connections} connections, but aria2 supports at most "
            f"{ARIA2_MAX_CONNECTIONS_PER_SERVER}.  Using {ARIA2_MAX_CONNECTIONS_PER_SERVER}."
        )
        return ARIA2_MAX_CONNECTIONS_PER_SERVER
    return requested_connections
