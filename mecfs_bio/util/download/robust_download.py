import tempfile
import time
from pathlib import Path
from subprocess import CalledProcessError

import structlog

from mecfs_bio.util.download.verify import hash_matches
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()


def robust_download_with_curl(
    md5sum: str | None,
    dest: Path,
    url: str,
    max_outer_retries: int = 10,
):
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        temp_out = tmp_path / dest.name
        for i in range(max_outer_retries):
            try:
                execute_command(
                    cmd=[
                        "curl",
                        "-L",
                        "-C",
                        "-",
                        "--retry",
                        "8",
                        "--retry-all-errors",
                        "--retry-delay",
                        "5",
                        "--connect-timeout",
                        "10",
                        "--max-time",
                        "0",
                        "-o",
                        str(temp_out),
                        url,
                    ]
                )
                if temp_out.exists() and hash_matches(temp_out, md5sum):
                    temp_out.rename(dest)
                    return

            except CalledProcessError as e:
                if i >= (max_outer_retries - 1):
                    break
                backoff = min(2 ** (i), 60)
                logger.debug(
                    f"Download attempt {i + 1} failed.  Backing off for {backoff} seconds"
                )
                time.sleep(backoff)
        raise RuntimeError(
            f"Download of from {url} to {dest} failed after {max_outer_retries}"
        )
