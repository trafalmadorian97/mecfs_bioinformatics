"""
Pytest plugin for measuring the cost of explicit gc.collect() calls during the unit-test run.

Wraps gc.collect so every call (from gwaslab, from mecfs_bio tasks and the runner's post-execute
hook, from pytest) is timed, and prints the total at the end of the session. With GC_FREEZE=1 it
also calls gc.freeze() once collection finishes, which moves every object alive at that point
(modules, test items, fixtures definitions) into a permanent generation that later full
collections skip.

Experiment-only: load with PYTHONPATH=experiments/claude/test_speed pytest -p gc_timing_plugin.
"""

import gc
import os
import time

_original_collect = gc.collect
_stats = {"calls": 0, "seconds": 0.0}


def _timed_collect(*args: int) -> int:
    start = time.perf_counter()
    result = _original_collect(*args)
    _stats["seconds"] += time.perf_counter() - start
    _stats["calls"] += 1
    return result


gc.collect = _timed_collect


def pytest_collection_finish(session: object) -> None:
    if os.environ.get("GC_FREEZE") == "1":
        gc.freeze()


def pytest_sessionfinish(session: object) -> None:
    worker = os.environ.get("PYTEST_XDIST_WORKER", "main")
    print(
        f"\nGCSTATS worker={worker} freeze={os.environ.get('GC_FREEZE', '0')} "
        f"calls={_stats['calls']} seconds={_stats['seconds']:.2f}",
        flush=True,
    )
