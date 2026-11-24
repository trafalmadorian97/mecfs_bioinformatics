from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.wf.parallel_range_downloader import (
    HttpParallelRangeDownloader,
    ParallelRangeDownloader,
)
from mecfs_bio.build_system.wf.s3_access import Boto3S3Access, S3Access
from mecfs_bio.build_system.wf.synapse_downloader import (
    SharedClientSynapseDownloader,
    SynapseDownloader,
)
from mecfs_bio.build_system.wf.wf_downloader import (
    SimpleWFDownloader,
    WFDownloader,
)
from mecfs_bio.util.download.parallel_ranges import MIB
from mecfs_bio.util.file_io.ranged_tar import RangeReader


@frozen(slots=True)
class WF:
    """
    An interface to the external world.

    To add a new capability, add an attribute here and a corresponding parameter to make_wf
    with an appropriate default.
    """

    downloader: WFDownloader
    synapse_downloader: SynapseDownloader
    s3_access: S3Access
    parallel_range_downloader: ParallelRangeDownloader

    def download_from_url(
        self,
        url: str,
        local_path: Path,
        md5_hash: str | None,
        requested_connections: int | None = None,
    ) -> None:
        """
        Download the file at url into local_path and verify its md5 hash.

        requested_connections is a hint for the number of simultaneous connections to use.
        Not all downloaders honor it.
        """
        self.downloader.download(
            url=url,
            local_path=local_path,
            md5_hash=md5_hash,
            requested_connections=requested_connections,
        )

    def download_url_parallel(
        self,
        url: str,
        size: int,
        local_path: Path,
        md5_hash: str,
        n_connections: int = 128,
        chunk_bytes: int = 4 * MIB,
    ) -> None:
        """
        Download the size bytes at url into local_path with n_connections parallel HTTP Range
        requests of chunk_bytes each, then verify the md5 hash.
        """
        self.parallel_range_downloader.download(
            url=url,
            size=size,
            local_path=local_path,
            md5_hash=md5_hash,
            n_connections=n_connections,
            chunk_bytes=chunk_bytes,
        )

    def download_from_synapse(
        self, synid: str, dest_dir: Path, expected_name: str
    ) -> Path:
        """Download the single file at synid into dest_dir and return its path."""
        return self.synapse_downloader.download(synid, dest_dir, expected_name)

    def fetch_synapse_file_head(self, synid: str, n_bytes: int) -> bytes:
        """Return the first n_bytes of the Synapse file at synid without downloading the
        whole file (via a pre-signed URL + HTTP Range request)."""
        return self.synapse_downloader.fetch_file_head(synid, n_bytes)

    def open_synapse_range_reader(self, synid: str) -> RangeReader:
        """A reader of byte ranges of the Synapse file at synid, for a burst of reads that
        avoids downloading the whole file."""
        return self.synapse_downloader.open_range_reader(synid)


def make_wf(
    downloader: WFDownloader | None = None,
    synapse_downloader: SynapseDownloader | None = None,
    s3_access: S3Access | None = None,
    parallel_range_downloader: ParallelRangeDownloader | None = None,
) -> WF:
    """Construct a WF, defaulting any capability not explicitly supplied."""
    return WF(
        downloader=downloader if downloader is not None else SimpleWFDownloader(),
        synapse_downloader=(
            synapse_downloader
            if synapse_downloader is not None
            else SharedClientSynapseDownloader()
        ),
        s3_access=s3_access if s3_access is not None else Boto3S3Access(),
        parallel_range_downloader=(
            parallel_range_downloader
            if parallel_range_downloader is not None
            else HttpParallelRangeDownloader()
        ),
    )
