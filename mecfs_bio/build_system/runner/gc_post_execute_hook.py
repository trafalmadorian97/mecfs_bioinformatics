"""
Post-execute hook that runs Python's cyclic garbage collector between tasks.

Long-running gwaslab pipelines accumulate cycle-held objects (Sumstats/LDSC
internals) that refcounting alone cannot reclaim. Forcing a `gc.collect()`
after each task reclaims multi-GB of resident memory per cell-analysis task
in the s-LDSC system test.
"""

import gc

import structlog

from mecfs_bio.build_system.task.base_task import Task

logger = structlog.get_logger()


def _current_rss_kb() -> int:
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except OSError:
        pass
    return -1


def gc_collect_post_execute_hook(task: Task) -> None:
    rss_before_kb = _current_rss_kb()
    collected = gc.collect()
    rss_after_kb = _current_rss_kb()
    logger.debug(
        f"post-task gc.collect for {task.asset_id}: "
        f"unreachable_collected={collected} "
        f"rss_before_kb={rss_before_kb} rss_after_kb={rss_after_kb}"
    )
