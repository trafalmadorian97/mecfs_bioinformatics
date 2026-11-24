"""
WF capability: download a large file with many parallel HTTP Range requests, then verify md5.
"""

from abc import ABC, abstractmethod
from pathlib import Path

from attrs import frozen

from mecfs_bio.util.download.parallel_ranges import (
    REQUESTS_RANGE_GETTER,
    RangeGetter,
    download_with_parallel_ranges,
)
from mecfs_bio.util.download.verify import verify_hash


class ParallelRangeDownloader(ABC):
    @abstractmethod
    def download(
        self,
        url: str,
        size: int,
        local_path: Path,
        md5_hash: str,
        n_connections: int,
        chunk_bytes: int,
    ) -> None:
        """Download the size bytes of url into local_path and verify their md5."""


@frozen(slots=True)
class HttpParallelRangeDownloader(ParallelRangeDownloader):
    getter: RangeGetter = REQUESTS_RANGE_GETTER

    def download(
        self,
        url: str,
        size: int,
        local_path: Path,
        md5_hash: str,
        n_connections: int,
        chunk_bytes: int,
    ) -> None:
        download_with_parallel_ranges(
            url,
            size,
            local_path,
            n_connections=n_connections,
            chunk_bytes=chunk_bytes,
            getter=self.getter,
        )
        verify_hash(local_path, md5_hash)
