from abc import ABC
from pathlib import Path

import py3_wget
from attrs import Factory, frozen

from mecfs_bio.build_system.wf.synapse_downloader import (
    SharedClientSynapseDownloader,
    SynapseDownloader,
)
from mecfs_bio.util.download.robust_download import robust_download_with_aria
from mecfs_bio.util.download.verify import verify_hash


@frozen
class WF(ABC):
    """
    An interface to the external world.
    """

    synapse_downloader: SynapseDownloader = Factory(SharedClientSynapseDownloader)

    def download_from_url(
        self, url: str, local_path: Path, md5_hash: str | None
    ) -> None:
        pass

    def download_from_synapse(
        self, synid: str, dest_dir: Path, expected_name: str
    ) -> Path:
        """Download the single file at synid into dest_dir and return its path."""
        return self.synapse_downloader.download(synid, dest_dir, expected_name)


@frozen
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


@frozen
class RobustDownloadWF(WF):
    def download_from_url(
        self, url: str, local_path: Path, md5_hash: str | None
    ) -> None:
        robust_download_with_aria(
            md5sum=md5_hash,
            url=url,
            dest=local_path,
        )
