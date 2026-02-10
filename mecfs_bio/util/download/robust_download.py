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
    num_simil: int=4,
):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        temp_out = tmp_path / dest.name
        last_attempt_succeeded = False
        for i in range(max_outer_retries):
            if (
                temp_out.exists()
                and hash_matches(temp_out, md5sum)
                and last_attempt_succeeded
            ):
                temp_out.rename(dest)
                return
            try:
                cmd = [
                    "stdbuf", "-o0", "-e0",
                    "pixi",
                    "r",
                    "--environment",
                    "download-env",
                    "aria2c",
                    "-x",
                    str(num_simil),
                    "--continue=true",
                    "--check-integrity=true",
                    "--allow-overwrite=true",
                    "--auto-file-renaming=false",
                    "--max-tries=8",
                    "--retry-wait=5",
                    "--timeout=30",
                    "--connect-timeout=10",
                    "--file-allocation=none",
                    "--dir", str(temp_out.parent),
                    "--out", temp_out.name,
                    url,
                ]
                # --- DEFINE ENVIRONMENTS ---
                execute_command(
                    cmd=cmd
                )
                last_attempt_succeeded = True
            except CalledProcessError as e:
                if i >= max_outer_retries:
                    break
                last_attempt_succeeded = False
                backoff = min(2 ** (i), 60)
                logger.debug(
                    f"Download attempt {i + 1} failed.  Backing off for {backoff} seconds"
                )
                time.sleep(backoff)
        raise RuntimeError(
            f"Download of from {url} to {dest} failed after {max_outer_retries}"
        )

