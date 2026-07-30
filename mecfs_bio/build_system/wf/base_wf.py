from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.wf.synapse_downloader import (
    SharedClientSynapseDownloader,
    SynapseDownloader,
)
from mecfs_bio.build_system.wf.wf_downloader import (
    SimpleWFDownloader,
    WFDownloader,
)


@frozen
class WF:
    """
    An interface to the external world.

    To add a new capability, add an attribute here and a corresponding parameter to make_wf
    with an appropriate default.
    """

    downloader: WFDownloader
    synapse_downloader: SynapseDownloader

    def download_from_url(
        self, url: str, local_path: Path, md5_hash: str | None
    ) -> None:
        self.downloader.download(url, local_path, md5_hash)

    def download_from_synapse(
        self, synid: str, dest_dir: Path, expected_name: str
    ) -> Path:
        """Download the single file at synid into dest_dir and return its path."""
        return self.synapse_downloader.download(synid, dest_dir, expected_name)

    def fetch_synapse_file_head(self, synid: str, n_bytes: int) -> bytes:
        """Return the first n_bytes of the Synapse file at synid without downloading the
        whole file (via a pre-signed URL + HTTP Range request)."""
        return self.synapse_downloader.fetch_file_head(synid, n_bytes)


def make_wf(
    downloader: WFDownloader | None = None,
    synapse_downloader: SynapseDownloader | None = None,
) -> WF:
    """Construct a WF, defaulting any capability not explicitly supplied."""
    return WF(
        downloader=downloader if downloader is not None else SimpleWFDownloader(),
        synapse_downloader=(
            synapse_downloader
            if synapse_downloader is not None
            else SharedClientSynapseDownloader()
        ),
    )
