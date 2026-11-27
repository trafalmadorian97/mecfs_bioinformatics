import hashlib
import random
import threading
from pathlib import Path

import pytest
from attrs import define, field

from mecfs_bio.build_system.wf.parallel_range_downloader import (
    HttpParallelRangeDownloader,
)
from mecfs_bio.util.download.parallel_ranges import (
    MIB,
    RangeGetter,
    RangeResponse,
    download_with_parallel_ranges,
)

_SOURCE = random.Random(0).randbytes(10 * MIB)


def _no_sleep(seconds: float) -> None:
    pass


@define(slots=True)
class _BytesGetter(RangeGetter):
    """Serves ranges of _SOURCE; fails_first_per_chunk returns 503 on each chunk's first try."""

    status: int = 206
    fails_first_per_chunk: bool = False
    _seen: set[int] = field(factory=set)
    _lock: threading.Lock = field(factory=threading.Lock)

    def get(self, url: str, start: int, end_inclusive: int) -> RangeResponse:
        if self.fails_first_per_chunk:
            with self._lock:
                first = start not in self._seen
                self._seen.add(start)
            if first:
                return RangeResponse(status=503, body=b"")
        return RangeResponse(
            status=self.status, body=_SOURCE[start : end_inclusive + 1]
        )


def test_reassembles_exact_bytes(tmp_path: Path) -> None:
    out = tmp_path / "out.bin"
    download_with_parallel_ranges(
        "https://x",
        len(_SOURCE),
        out,
        n_connections=8,
        chunk_bytes=MIB + 7,
        getter=_BytesGetter(),
    )
    assert out.read_bytes() == _SOURCE


def test_retries_transient_5xx(tmp_path: Path) -> None:
    out = tmp_path / "out.bin"
    download_with_parallel_ranges(
        "https://x",
        len(_SOURCE),
        out,
        n_connections=8,
        chunk_bytes=MIB + 7,
        getter=_BytesGetter(fails_first_per_chunk=True),
        retry_sleep=_no_sleep,
    )
    assert out.read_bytes() == _SOURCE


def test_rejects_200_full_body(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        download_with_parallel_ranges(
            "https://x",
            len(_SOURCE),
            tmp_path / "out.bin",
            n_connections=2,
            chunk_bytes=4 * MIB,
            getter=_BytesGetter(status=200),
        )


def test_md5_mismatch_fails(tmp_path: Path) -> None:
    downloader = HttpParallelRangeDownloader(getter=_BytesGetter())
    out = tmp_path / "out.bin"
    good = hashlib.md5(_SOURCE).hexdigest()
    downloader.download(
        "https://x",
        len(_SOURCE),
        out,
        md5_hash=good,
        n_connections=4,
        chunk_bytes=2 * MIB,
    )
    with pytest.raises(AssertionError):
        downloader.download(
            "https://x",
            len(_SOURCE),
            out,
            md5_hash="0" * 32,
            n_connections=4,
            chunk_bytes=2 * MIB,
        )
