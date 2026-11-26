import os
import subprocess
import sys
from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.util.subproc.run_in_child_process import run_in_child_process


def _pid_and_sum(a: int, b: int) -> tuple[int, int]:
    return os.getpid(), a + b


def _fail(message: str) -> None:
    raise ValueError(message)


def _double_and_sum(frame: pl.DataFrame) -> pl.DataFrame:
    return frame.with_columns(pl.col("x") * 2).select(pl.col("x").sum())


def _exit_hard() -> None:
    os._exit(3)


def test_runs_in_another_process_and_returns_result() -> None:
    pid, total = run_in_child_process(_pid_and_sum, 2, b=3)
    assert pid != os.getpid()
    assert total == 5


def test_reraises_child_exception() -> None:
    with pytest.raises(ValueError):
        run_in_child_process(_fail, "boom")


def test_child_runs_polars_while_parent_has_polars_threads() -> None:
    """The parent has run polars (so its thread pool exists); the child must not hang."""
    frame = pl.DataFrame({"x": range(100_000)})
    assert frame.select(pl.col("x").sum()).item() > 0
    result = run_in_child_process(_double_and_sum, frame)
    assert result.item() == 2 * sum(range(100_000))


_POLARS_MAIN = """
import polars as pl
from mecfs_bio.util.subproc.run_in_child_process import run_in_child_process

# Module-level work, as asset modules do when imported (e.g. reading a manifest).
frame = pl.DataFrame({"x": range(1_000_000)}).with_columns(y=pl.col("x") * 2)
frame.group_by(pl.col("x") % 7).agg(pl.col("y").sum())

if __name__ == "__main__":
    print(run_in_child_process(pl.DataFrame.sum, frame).item(0, 1))
"""


def test_child_does_not_hang_when_main_module_uses_polars(tmp_path: Path) -> None:
    """
    A build's main module runs polars when imported, which starts its thread pool. A child
    forked from a process holding those threads can deadlock; this one must finish.
    """
    script = tmp_path / "polars_main.py"
    script.write_text(_POLARS_MAIN)
    done = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
        timeout=30,
        check=True,
        cwd=Path(__file__).parents[4],
    )
    assert int(done.stdout.split()[-1]) == 2 * sum(range(1_000_000))


def test_child_that_dies_raises() -> None:
    with pytest.raises(ChildProcessError):
        run_in_child_process(_exit_hard)
