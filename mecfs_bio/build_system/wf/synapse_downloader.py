"""
Synapse download capability
"""

from __future__ import annotations

import threading
from abc import ABC, abstractmethod
from functools import partial
from pathlib import Path

import requests
import synapseclient
import synapseutils
from attrs import frozen

from mecfs_bio.util.file_io.ranged_tar import RangeReader
from mecfs_bio.util.retry.call_with_retries import call_with_retries

_HTTP_PARTIAL_CONTENT = 206
_RANGE_TIMEOUT_SECONDS = 120


class SynapseDownloader(ABC):
    """Downloads a single file from synapse.org into a local directory."""

    @abstractmethod
    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        """Download the single file at synid into dest_dir, assert its name equals
        expected_name, and return its local path."""

    @abstractmethod
    def open_range_reader(self, synid: str) -> RangeReader:
        """A reader of arbitrary byte ranges of the file at synid, which never downloads
        the whole file. Meant for a burst of reads soon after opening: the underlying link
        may expire."""

    def fetch_file_head(self, synid: str, n_bytes: int) -> bytes:
        """Return the first n_bytes of the file at synid WITHOUT downloading the whole
        file. Used to read a tiny prefix (e.g. a tar header + the first rows of the first
        gzip member) from a large entity."""
        return self.open_range_reader(synid).read_range(0, n_bytes)


class _RetryableHttpError(Exception):
    """A server-side (5xx) HTTP failure, which a later identical request may survive."""


class _ThreadLocalSessions:
    """
    One requests.Session per thread, created on first use.

    A session keeps its HTTPS connections open between requests, so consecutive range reads
    of the same host skip connection setup (several network round trips each). The
    pre-signed URLs of all Synapse files point at the same S3 host, so the connections also
    carry over from one file to the next. requests.Session is not guaranteed thread-safe,
    hence one per thread.
    """

    def __init__(self) -> None:
        self._local = threading.local()

    def get(self) -> requests.Session:
        session: requests.Session | None = getattr(self._local, "session", None)
        if session is None:
            session = requests.Session()
            self._local.session = session
        return session


@frozen(slots=True)
class _PresignedUrlRangeReader:
    """Reads byte ranges of a pre-signed S3 URL with HTTP Range requests, each on the
    calling thread's session."""

    presigned_url: str
    sessions: _ThreadLocalSessions

    def _get(self, start: int, length: int) -> bytes:
        response = self.sessions.get().get(
            self.presigned_url,
            headers={"Range": f"bytes={start}-{start + length - 1}"},
            timeout=_RANGE_TIMEOUT_SECONDS,
        )
        if response.status_code >= 500:
            raise _RetryableHttpError(f"HTTP {response.status_code} for a range read")
        response.raise_for_status()
        # A 200 would mean the server ignored Range and sent the whole file.
        assert response.status_code == _HTTP_PARTIAL_CONTENT, (
            f"expected HTTP 206 for a range read, got {response.status_code}"
        )
        return response.content

    def read_range(self, start: int, length: int) -> bytes:
        assert start >= 0 and length >= 1, f"bad range: start={start} length={length}"
        return call_with_retries(
            partial(self._get, start, length),
            retry_on=(
                _RetryableHttpError,
                requests.ConnectionError,
                requests.Timeout,
                requests.exceptions.ChunkedEncodingError,
            ),
        )


class SharedClientSynapseDownloader(SynapseDownloader):
    """Downloads via a lazily authenticated Synapse client reused across calls."""

    def __init__(self) -> None:
        self._client: synapseclient.Synapse | None = None
        self._login_lock = threading.Lock()
        self._sessions = _ThreadLocalSessions()

    def _client_or_login(self) -> synapseclient.Synapse:
        with self._login_lock:
            if self._client is None:
                self._client = synapseclient.login()
            return self._client

    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:

        files = synapseutils.syncFromSynapse(
            self._client_or_login(), synid, path=str(dest_dir)
        )
        assert len(files) == 1, f"expected one file for {synid}, got {len(files)}"
        assert files[0].name == expected_name, files[0].name
        return Path(files[0].path)

    def open_range_reader(self, synid: str) -> RangeReader:
        # Synapse entities are S3-backed; the file handle yields a pre-signed S3 URL that
        # honours HTTP Range, so parts of an (often several-hundred-MB) file can be read
        # without downloading all of it.
        client = self._client_or_login()
        entity = client.get(synid, downloadFile=False)
        download_info = client._getFileHandleDownload(entity.dataFileHandleId, synid)
        return _PresignedUrlRangeReader(
            presigned_url=download_info["preSignedURL"], sessions=self._sessions
        )
