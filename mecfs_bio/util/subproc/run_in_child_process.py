"""
Run one function call in a short-lived child process, so that memory the call leaks (for
example inside a native library) is returned to the system when the child exits.

The child is a fresh Python interpreter (python -m child_call_main), not a fork. Forking a
process that runs polars can deadlock the child: polars, its jemalloc allocator and its tokio
runtime start threads, and a fork copies their locks but not the threads that would release
them. A fresh interpreter also avoids multiprocessing's spawn, which re-imports the parent's
main module in the child (seconds for a build entry point that imports the asset graph).

The function, its arguments and its result are passed as pickle files in a temporary
directory, so all three must be picklable and the function must be importable by name (a
module-level function). An exception raised in the child is re-raised in the parent; a child
that dies without writing a result raises ChildProcessError.
"""

import pickle
import subprocess
import sys
import tempfile
from collections.abc import Callable
from pathlib import Path
from typing import cast

from mecfs_bio.util.subproc import child_call, child_call_main


def run_in_child_process[**P, R](
    fn: Callable[P, R], *args: P.args, **kwargs: P.kwargs
) -> R:
    with tempfile.TemporaryDirectory() as tmp:
        call_path = Path(tmp) / child_call.CALL_FILENAME
        result_path = Path(tmp) / child_call.RESULT_FILENAME
        with call_path.open("wb") as f:
            pickle.dump((fn, args, kwargs), f, protocol=pickle.HIGHEST_PROTOCOL)
        completed = subprocess.run(
            [sys.executable, "-m", child_call_main.__name__, tmp], check=False
        )
        if not result_path.exists():
            raise ChildProcessError(
                f"child process for {fn!r} exited with code "
                f"{completed.returncode} without a result"
            )
        with result_path.open("rb") as f:
            outcome = pickle.load(f)
    if isinstance(outcome, child_call.ChildFailure):
        raise outcome.exception
    assert isinstance(outcome, child_call.ChildSuccess)
    return cast(R, outcome.value)
