"""
The child side of run_in_child_process: run_call reads the pickled (function, args, kwargs)
from a directory and writes the outcome to the same directory: ChildSuccess with the return
value, or ChildFailure with the exception. An exception that cannot be pickled is replaced by a
RuntimeError carrying its formatted traceback. The child's entry point is child_call_main.
"""

import pickle
import traceback
from pathlib import Path

from attrs import frozen

CALL_FILENAME = "call.pkl"
RESULT_FILENAME = "result.pkl"


@frozen(slots=True)
class ChildSuccess:
    value: object


@frozen(slots=True)
class ChildFailure:
    exception: BaseException


def _write(outcome: ChildSuccess | ChildFailure, path: Path) -> None:
    with path.open("wb") as f:
        pickle.dump(outcome, f, protocol=pickle.HIGHEST_PROTOCOL)


def run_call(call_dir: Path) -> None:
    with (call_dir / CALL_FILENAME).open("rb") as f:
        fn, args, kwargs = pickle.load(f)
    result_path = call_dir / RESULT_FILENAME
    try:
        outcome: ChildSuccess | ChildFailure = ChildSuccess(fn(*args, **kwargs))
    except Exception as e:
        outcome = ChildFailure(e)
    try:
        _write(outcome, result_path)
    except (pickle.PicklingError, TypeError, AttributeError):
        assert isinstance(outcome, ChildFailure)
        formatted = "".join(traceback.format_exception(outcome.exception))
        _write(ChildFailure(RuntimeError(formatted)), result_path)
