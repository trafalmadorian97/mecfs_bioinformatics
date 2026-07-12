"""
Synapse download capability for the WF external-world interface.

WF delegates download_from_synapse to a SynapseDownloader so that the delivery
strategy (how a client is obtained and reused) is composed into WF rather than
baked into each WF subclass. The default, SharedClientSynapseDownloader, logs in
once and reuses the client across a run so a pipeline that pulls many Synapse
entities does not re-authenticate per task.
"""

from __future__ import annotations

import threading
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import synapseclient


class SynapseDownloader(ABC):
    """Downloads a single file from synapse.org into a local directory."""

    @abstractmethod
    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        """Download the single file at synid into dest_dir, assert its name equals
        expected_name, and return its local path."""


class SharedClientSynapseDownloader(SynapseDownloader):
    """Downloads via a lazily authenticated Synapse client reused across calls.

    The client is created on first use behind a lock, so a run of many downloads
    authenticates once. The lock guards lazy login only; relying on a single
    client for concurrent syncFromSynapse calls is deferred until the scheduler
    actually runs tasks in parallel.
    """

    def __init__(self) -> None:
        self._client: synapseclient.Synapse | None = None
        self._login_lock = threading.Lock()

    def _client_or_login(self) -> synapseclient.Synapse:
        with self._login_lock:
            if self._client is None:
                import synapseclient

                self._client = synapseclient.login()
            return self._client

    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        import synapseutils

        files = synapseutils.syncFromSynapse(
            self._client_or_login(), synid, path=str(dest_dir)
        )
        assert len(files) == 1, f"expected one file for {synid}, got {len(files)}"
        assert files[0].name == expected_name, files[0].name
        return Path(files[0].path)
