"""
Download one large file over HTTP with many parallel Range requests.

Some object stores serve each connection slowly (deCODE's: about 0.2 MB/s per connection) but
scale almost linearly with the number of connections. This downloader splits the file into
fixed-size chunks, fetches them from a thread pool, and writes each chunk at its offset in a
preallocated file. Each chunk is retried on transient failures (connection errors, timeouts,
truncated bodies, 5xx). A server that answers 200 with the whole body instead of 206 fails
the download rather than silently writing too much.
"""

import os
import threading
import time
from abc import ABC, abstractmethod
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path

import requests
from attrs import frozen
from requests.exceptions import ChunkedEncodingError

from mecfs_bio.util.retry.call_with_retries import call_with_retries

MIB = 1024 * 1024
_PARTIAL_CONTENT = 206
_TIMEOUT_S = 120


@frozen(slots=True)
class RangeResponse:
    status: int
    body: bytes


class RangeGetter(ABC):
    @abstractmethod
    def get(self, url: str, start: int, end_inclusive: int) -> RangeResponse:
        """GET bytes start..end_inclusive of url."""


class _ThreadSessions(threading.local):
    session: requests.Session | None = None


_SESSIONS = _ThreadSessions()


@frozen(slots=True)
class RequestsRangeGetter(RangeGetter):
    """Range GETs through one requests.Session per thread."""

    def get(self, url: str, start: int, end_inclusive: int) -> RangeResponse:
        if _SESSIONS.session is None:
            _SESSIONS.session = requests.Session()
        r = _SESSIONS.session.get(
            url, headers={"Range": f"bytes={start}-{end_inclusive}"}, timeout=_TIMEOUT_S
        )
        return RangeResponse(status=r.status_code, body=r.content)


REQUESTS_RANGE_GETTER = RequestsRangeGetter()


class _RetryableRangeError(Exception):
    pass


_RETRY_ON: tuple[type[Exception], ...] = (
    requests.ConnectionError,
    requests.Timeout,
    ChunkedEncodingError,
    _RetryableRangeError,
)


@frozen(slots=True)
class _Chunk:
    start: int
    end_inclusive: int

    @property
    def length(self) -> int:
        return self.end_inclusive - self.start + 1


def chunk_bounds(start: int, stop: int, chunk_bytes: int) -> list[_Chunk]:
    """Inclusive chunks covering [start, stop)."""
    return [
        _Chunk(s, min(s + chunk_bytes, stop) - 1)
        for s in range(start, stop, chunk_bytes)
    ]


def fetch_range(getter: RangeGetter, url: str, start: int, end_inclusive: int) -> bytes:
    """One range GET, raising a retryable error on a 5xx and asserting a 206 of exact length."""
    response = getter.get(url, start, end_inclusive)
    if response.status >= 500:
        raise _RetryableRangeError(
            f"HTTP {response.status} for bytes {start}-{end_inclusive}"
        )
    assert response.status == _PARTIAL_CONTENT, (
        f"expected HTTP 206 for bytes {start}-{end_inclusive}, got {response.status}"
    )
    expected = end_inclusive - start + 1
    if len(response.body) != expected:
        raise _RetryableRangeError(
            f"short body for bytes {start}-{end_inclusive}: {len(response.body)} of {expected}"
        )
    return response.body


def fetch_range_with_retries(
    getter: RangeGetter,
    url: str,
    start: int,
    end_inclusive: int,
    sleep: Callable[[float], None] = time.sleep,
) -> bytes:
    return call_with_retries(
        partial(fetch_range, getter, url, start, end_inclusive),
        retry_on=_RETRY_ON,
        max_attempts=8,
        sleep=sleep,
    )


def download_with_parallel_ranges(
    url: str,
    size: int,
    out_path: Path,
    n_connections: int = 128,
    chunk_bytes: int = 4 * MIB,
    getter: RangeGetter = REQUESTS_RANGE_GETTER,
    retry_sleep: Callable[[float], None] = time.sleep,
) -> None:
    """Download the size bytes of url into out_path, n_connections chunks at a time."""
    assert size >= 0 and n_connections >= 1 and chunk_bytes >= 1
    fd = os.open(out_path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    try:
        os.truncate(fd, size)

        def one(chunk: _Chunk) -> None:
            body = fetch_range_with_retries(
                getter, url, chunk.start, chunk.end_inclusive, sleep=retry_sleep
            )
            os.pwrite(fd, body, chunk.start)

        with ThreadPoolExecutor(n_connections) as pool:
            list(pool.map(one, chunk_bounds(0, size, chunk_bytes)))
    finally:
        os.close(fd)
