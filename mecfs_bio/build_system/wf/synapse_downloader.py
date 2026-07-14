"""
Synapse download capability
"""

from __future__ import annotations

import threading
from abc import ABC, abstractmethod
from pathlib import Path

import synapseclient
import synapseutils


class SynapseDownloader(ABC):
    """Downloads a single file from synapse.org into a local directory."""

    @abstractmethod
    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        """Download the single file at synid into dest_dir, assert its name equals
        expected_name, and return its local path."""


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
