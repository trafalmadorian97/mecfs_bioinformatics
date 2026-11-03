import tempfile
import time
from abc import ABC, abstractmethod
from collections.abc import Callable
from functools import partial
from pathlib import Path
from subprocess import CalledProcessError

import structlog
from attrs import frozen

from mecfs_bio.util.download.verify import hash_matches
from mecfs_bio.util.retry.call_with_retries import call_with_retries
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

ARIA2_MAX_CONNECTIONS_PER_SERVER = 16
_UNVERIFIED_FILE_PREVIEW_BYTES = 512


class DownloadAttemptFailedError(RuntimeError):
    pass


class Downloader(ABC):
    @abstractmethod
    def download(
        self,
        url: str,
        local_path: Path,
    ) -> bool:
        pass


def build_aria2c_command(
    url: str, local_path: Path, summary_interval: int, num_connections: int
) -> list[str]:
    """
    Build the aria2c command line that downloads url to local_path.

    num_connections is passed as both -x (max connections per server) and -s (split count),
    since aria2 opens min(-x, -s) connections to a single server.
    """
    return [
        "pixi",
        "r",
        "--environment",
        "download-env",
        "aria2c",
        f"--summary-interval={summary_interval}",
        "-x",
        str(num_connections),
        "-s",
        str(num_connections),
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


@frozen(slots=True)
class AriaDownloader(Downloader):
    """
    Downloader that uses aria2
    https://aria2.github.io/manual/en/html/aria2c.html

    num_simil is the number of simultaneous connections to the server, between 1 and
    ARIA2_MAX_CONNECTIONS_PER_SERVER.
    """

    summary_interval: int = 10
    num_simil: int = 1
    executor: Callable[[list[str]], str] = execute_command

    def __attrs_post_init__(self) -> None:
        assert 1 <= self.num_simil <= ARIA2_MAX_CONNECTIONS_PER_SERVER, (
            f"num_simil must be between 1 and {ARIA2_MAX_CONNECTIONS_PER_SERVER}, got {self.num_simil}"
        )

    def download(self, url: str, local_path: Path) -> bool:
        try:
            self.executor(
                build_aria2c_command(
                    url=url,
                    local_path=local_path,
                    summary_interval=self.summary_interval,
                    num_connections=self.num_simil,
                )
            )
            return True
        except CalledProcessError as e:
            logger.error(f"Failed to download {url}: {e}")
            return False


def _log_unverified_file(path: Path) -> None:
    """
    Log the size and leading bytes of a rejected download, since its temp directory is
    deleted afterwards and it cannot be inspected later.
    """
    if not path.exists():
        logger.warning(f"Downloader reported success, but {path} does not exist.")
        return
    with open(path, "rb") as f:
        leading_bytes = f.read(_UNVERIFIED_FILE_PREVIEW_BYTES)
    logger.warning(
        f"Rejected downloaded file {path}: size {path.stat().st_size} bytes, "
        f"leading bytes {leading_bytes!r}"
    )


def _attempt_download(
    md5sum: str | None,
    dest: Path,
    url: str,
    downloader: Downloader,
    temp_out: Path,
) -> None:
    """
    Make one download attempt, moving a verified file to dest.

    A rejected file is deleted: aria2 with --continue=true treats an existing file of the
    full size as finished, so keeping it would make the next attempt return the same bytes.
    """
    logger.debug(f"Downloading from {url} to {temp_out}")
    if not downloader.download(url, local_path=temp_out):
        raise DownloadAttemptFailedError(f"Downloader failed to fetch {url}.")
    if temp_out.exists() and hash_matches(temp_out, md5sum):
        temp_out.rename(dest)
        return
    _log_unverified_file(temp_out)
    temp_out.unlink(missing_ok=True)
    raise DownloadAttemptFailedError(f"Could not verify file downloaded from {url}.")


def robust_download(
    md5sum: str | None,
    dest: Path,
    url: str,
    downloader: Downloader,
    max_outer_retries: int = 10,
    max_backoff_time: int = 60,
    sleep: Callable[[float], None] = time.sleep,
) -> None:
    """
    Call a downloader repeatedly, with backoff, until it produces a verified file.

    A failed download and a file failing verification are retried identically: both can
    come from a transient server problem, such as an error page served with status 200.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        call_with_retries(
            partial(
                _attempt_download,
                md5sum=md5sum,
                dest=dest,
                url=url,
                downloader=downloader,
                temp_out=Path(tmpdir) / dest.name,
            ),
            retry_on=(DownloadAttemptFailedError,),
            max_attempts=max_outer_retries,
            max_backoff_seconds=max_backoff_time,
            sleep=sleep,
        )


def robust_download_with_aria(
    md5sum: str | None,
    dest: Path,
    url: str,
    max_outer_retries: int = 10,
    num_simil: int = 1,
    summary_interval: int = 10,
) -> None:
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
    )
