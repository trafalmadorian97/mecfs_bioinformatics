from subprocess import CalledProcessError

import pytest

from mecfs_bio.util.subproc.run_command import execute_command_with_retries


def test_retries_then_succeeds():
    """A command that fails twice then succeeds is retried to success."""
    calls = {"n": 0}
    slept: list[float] = []

    def flaky_executor(cmd: list[str]) -> str:
        calls["n"] += 1
        if calls["n"] < 3:
            raise CalledProcessError(returncode=1, cmd=cmd, output="HTTP 500")
        return "ok"

    result = execute_command_with_retries(
        ["echo", "hi"],
        base_backoff_seconds=0.0,
        sleep=slept.append,
        executor=flaky_executor,
    )

    assert result == "ok"
    assert calls["n"] == 3
    assert len(slept) == 2


def test_reraises_after_exhausting_attempts():
    """After max_attempts failures the last error is re-raised."""

    def always_fail(cmd: list[str]) -> str:
        raise CalledProcessError(returncode=1, cmd=cmd, output="HTTP 500")

    with pytest.raises(CalledProcessError):
        execute_command_with_retries(
            ["echo", "hi"],
            max_attempts=3,
            base_backoff_seconds=0.0,
            sleep=lambda _: None,
            executor=always_fail,
        )
