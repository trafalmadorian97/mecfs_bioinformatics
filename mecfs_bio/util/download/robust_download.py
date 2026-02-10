import tempfile
import time
from pathlib import Path
from subprocess import CalledProcessError

import structlog

from mecfs_bio.util.download.verify import hash_matches
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()


def robust_download_with_aria(
    md5sum: str | None,
    dest: Path,
    url: str,
    max_outer_retries: int = 10,
    num_simil: int = 4,
    summary_interval: int = 10,
):
    """
    Use aria2 to robustly download file.
    If aria2 fails, call it again in a loop
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        temp_out = tmp_path / dest.name
        for i in range(max_outer_retries):
            try:
                cmd = [
                    # "stdbuf", "-o0", "-e0",
                    "pixi",
                    "r",
                    "--environment",
                    "download-env",
                    "aria2c",
                    f"--summary-interval={summary_interval}",
                    "-x",
                    str(num_simil),
                    "--continue=true",
                    "--check-integrity=true",
                    "--allow-overwrite=true",
                    "--auto-file-renaming=false",
                    "--max-tries=8",
                    "--retry-wait=5",
                    "--timeout=30",
                    "--connect-timeout=30",
                    "--file-allocation=none",
                    "--dir",
                    str(temp_out.parent),
                    "--out",
                    temp_out.name,
                    url,
                ]
                execute_command(cmd=cmd)
                if temp_out.exists() and hash_matches(temp_out, md5sum):
                    temp_out.rename(dest)
                    return
            except CalledProcessError as e:
                if i >= (max_outer_retries-1):
                    break
                backoff = min(2 ** (i), 60)
                logger.debug(
                    f"Download attempt {i + 1} failed.  Backing off for {backoff} seconds"
                )
                time.sleep(backoff)
        raise RuntimeError(
            f"Download of from {url} to {dest} failed after {max_outer_retries}"
        )
