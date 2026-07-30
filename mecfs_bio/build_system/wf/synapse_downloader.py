"""
Synapse download capability
"""

from __future__ import annotations

import threading
from abc import ABC, abstractmethod
from pathlib import Path

import requests
import synapseclient
import synapseutils


class SynapseDownloader(ABC):
    """Downloads a single file from synapse.org into a local directory."""

    @abstractmethod
    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        """Download the single file at synid into dest_dir, assert its name equals
        expected_name, and return its local path."""

    @abstractmethod
    def fetch_file_head(self, synid: str, n_bytes: int) -> bytes:
        """Return the first n_bytes of the file at synid WITHOUT downloading the whole
        file. Used to read a tiny prefix (e.g. a tar header + the first rows of the first
        gzip member) from a large entity."""


class SharedClientSynapseDownloader(SynapseDownloader):
    """Downloads via a lazily authenticated Synapse client reused across calls."""

    def __init__(self) -> None:
        self._client: synapseclient.Synapse | None = None
        self._login_lock = threading.Lock()

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

    def fetch_file_head(self, synid: str, n_bytes: int) -> bytes:
        # Synapse entities are S3-backed; the file handle yields a pre-signed S3 URL that
        # honours HTTP Range, so we can pull just the first n_bytes instead of the whole
        # (often several-hundred-MB) file.
        client = self._client_or_login()
        entity = client.get(synid, downloadFile=False)
        file_handle_id = entity.dataFileHandleId
        download_info = client._getFileHandleDownload(file_handle_id, synid)
        presigned_url = download_info["preSignedURL"]
        response = requests.get(
            presigned_url,
            headers={"Range": f"bytes=0-{n_bytes - 1}"},
            timeout=60,
        )
        response.raise_for_status()
        return response.content
