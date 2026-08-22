import tempfile
import time
from abc import ABC, abstractmethod
from pathlib import Path
from subprocess import CalledProcessError

import structlog
from attrs import frozen

from mecfs_bio.util.download.verify import hash_matches
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()


class Downloader(ABC):
    @abstractmethod
    def download(
        self,
        url: str,
        local_path: Path,
        request_connections: int | None = None,
    ) -> bool:
        """Download url to local_path.

        request_connections is an advisory request for how many parallel connections
        to open for this one file, for callers that know the host tolerates it (some
        hosts throttle per connection, so more connections download a large file much
        faster). It is only a request: a downloader that cannot vary its connection
        count is free to ignore it. None means no preference.
        """


@frozen
class AriaDownloader(Downloader):
    """
    Downloader that uses aria2
    https://aria2.github.io/manual/en/html/aria2c.html
    """

    summary_interval: int = 10
    num_simil: int = 1

    def download(
        self, url: str, local_path: Path, request_connections: int | None = None
    ) -> bool:
        # Honor a per-file connection request; fall back to this downloader's default.
        connections = (
            request_connections if request_connections is not None else self.num_simil
        )
        cmd = aria_command(
            url=url,
            local_path=local_path,
            connections=connections,
            summary_interval=self.summary_interval,
        )
        try:
            execute_command(cmd=cmd)
            return True
        except CalledProcessError as e:
            logger.error(f"Failed to download {url}: {e}")
            return False


_ARIA_MAX_CONNECTIONS = 16


def aria_command(
    url: str, local_path: Path, connections: int, summary_interval: int
) -> list[str]:
    """Build the aria2c command line for downloading url to local_path.

    -s must match -x (and min-split-size be small enough) or aria2 caps the actual
    connection count regardless of -x. aria2 rejects --max-connection-per-server above
    16, so connections is asserted in range rather than left to fail mid-download.
    """
    assert 1 <= connections <= _ARIA_MAX_CONNECTIONS, (
        f"aria2 allows 1..{_ARIA_MAX_CONNECTIONS} connections per server, got {connections}"
    )
    return [
        "pixi",
        "r",
        "--environment",
        "download-env",
        "aria2c",
        f"--summary-interval={summary_interval}",
        "-x",
        str(connections),
        "-s",
        str(connections),
        "--min-split-size=1M",
        "--continue=true",
        "--allow-overwrite=true",
        "--user-agent=Wget/1.21.4",  # This is needed, otherwise Dropbox rejects download attempts
        "--auto-file-renaming=false",
        "--max-tries=8",
        "--retry-wait=5",
        "--timeout=30",
        "--connect-timeout=30",
        "--file-allocation=none",
        "--dir",
        str(local_path.parent),
        "--out",
        local_path.name,
        url,
    ]


def robust_download(
    md5sum: str | None,
    dest: Path,
    url: str,
    downloader: Downloader,
    max_outer_retries: int = 10,
    max_backoff_time: int = 60,
    request_connections: int | None = None,
):
    """
    Call a downloader in a loop, to add robustness
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        temp_out = tmp_path / dest.name
        for i in range(max_outer_retries):
            logger.debug(f"Downloading from {url} to {temp_out}")
            success = downloader.download(
                url, local_path=temp_out, request_connections=request_connections
            )
            if success:
                if temp_out.exists() and hash_matches(temp_out, md5sum):
                    temp_out.rename(dest)
                    return
                else:
                    logger.debug(
                        "Downloader returned success, but downloaded file could not be verified."
                    )
            else:
                if i >= (max_outer_retries - 1):
                    break
                backoff = min(2 ** (i), max_backoff_time)
                logger.debug(
                    f"Download attempt {i + 1} failed.  Backing off for {backoff} seconds."
                )
                time.sleep(backoff)
        raise RuntimeError(
            f"Download of from {url} to {dest} failed after {max_outer_retries} retries."
        )


def robust_download_with_aria(
    md5sum: str | None,
    dest: Path,
    url: str,
    max_outer_retries: int = 10,
    num_simil: int = 1,
    summary_interval: int = 10,
    request_connections: int | None = None,
):
    """
    Use aria2 to robustly download file.
    If aria2 fails, call it again in a loop

    https://aria2.github.io/manual/en/html/aria2c.html
    """
    robust_download(
        md5sum=md5sum,
        dest=dest,
        url=url,
        downloader=AriaDownloader(
            summary_interval=summary_interval, num_simil=num_simil
        ),
        max_outer_retries=max_outer_retries,
        request_connections=request_connections,
    )
