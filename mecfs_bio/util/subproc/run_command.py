import random
import subprocess
import sys
import time
from collections.abc import Callable
from subprocess import CalledProcessError

import structlog

logger = structlog.get_logger()


def execute_command(cmd: list[str]) -> str:
    """
    Execute and watch shell command, and return its output
    Source: https://stackoverflow.com/questions/29558074/output-of-subprocess-both-to-pipe-and-directly-to-stdout/77111541#77111541
    """
    for item in cmd:
        assert isinstance(item, str), f"{item} is not a string"
    output = ""
    str_cmd = " ".join(cmd)
    logger.debug(f"Running {str_cmd}")
    with subprocess.Popen(
        str_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        encoding="utf-8",
        errors="replace",
    ) as p:
        while True:
            stdout = p.stdout
            assert stdout is not None
            line = stdout.readline()
            if line != "":
                output += line + "\n"
                logger.debug(line)
            elif p.poll() != None:
                break
        sys.stdout.flush()
        if p.returncode == 0:
            return output
        else:
            raise CalledProcessError(cmd=cmd, returncode=p.returncode, output=output)


def execute_command_with_retries(
    cmd: list[str],
    max_attempts: int = 6,
    base_backoff_seconds: float = 2.0,
    max_backoff_seconds: float = 60.0,
    sleep: Callable[[float], None] = time.sleep,
    executor: Callable[[list[str]], str] = execute_command,
) -> str:
    """
    Run a command, retrying on CalledProcessError with exponential backoff
    and jitter.

    Intended for commands that hit flaky remote services (e.g. the GitHub
    release API, which intermittently returns HTTP 500 on asset downloads).
    The final failure is re-raised so callers still fail loudly on a genuinely
    broken command. The sleep and executor functions are injectable so tests
    can supply their own implementations without real delays or subprocesses.
    """
    assert max_attempts >= 1, f"max_attempts must be >= 1, got {max_attempts}"
    last_error: CalledProcessError | None = None
    for attempt in range(max_attempts):
        try:
            return executor(cmd)
        except CalledProcessError as e:
            last_error = e
            if attempt >= max_attempts - 1:
                break
            capped = min(base_backoff_seconds * 2**attempt, max_backoff_seconds)
            backoff = random.uniform(0, capped)
            logger.debug(
                f"Command failed (attempt {attempt + 1}/{max_attempts}, "
                f"exit {e.returncode}). Backing off {backoff:.1f}s before retry."
            )
            sleep(backoff)
    assert last_error is not None
    raise last_error
