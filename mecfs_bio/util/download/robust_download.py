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
    ) -> bool:
        pass


@frozen
class AriaDownloader(Downloader):
    """
    Downloader that uses aria2
    https://aria2.github.io/manual/en/html/aria2c.html
    """

    summary_interval: int = 10
    num_simil: int = 1

    def download(self, url: str, local_path: Path) -> bool:
        try:
            cmd = [
                "pixi",
                "r",
                "--environment",
                "download-env",
                "aria2c",
                f"--summary-interval={self.summary_interval}",
                "-x",
                str(self.num_simil),
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
            execute_command(cmd=cmd)
            return True
        except CalledProcessError as e:
            logger.error(f"Failed to download {url}: {e}")
            return False


def robust_download(
    md5sum: str | None,
    dest: Path,
    url: str,
    downloader: Downloader,
    max_outer_retries: int = 10,
    max_backoff_time: int = 60,
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
            success = downloader.download(url, local_path=temp_out)
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
    )
