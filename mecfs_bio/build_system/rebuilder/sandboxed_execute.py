import gc
import shutil
import tempfile
from pathlib import Path

import attr
import structlog
from loguru import logger


def _current_rss_kb() -> int:
    try:
        with open("/proc/self/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except OSError:
        pass
    return -1

logger = structlog.get_logger()

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.rebuilder.metadata_to_path.base_meta_to_path import (
    MetaToPath,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


def sandboxed_execute(
    task: Task,
    meta_to_path: MetaToPath,
    wf: WF,
    fetch: Fetch,
) -> Asset:
    """
    Execute a task in a temporary directory, then move its result to the final location
    determined by the task's metadata.
    """
    meta = task.meta
    target_path = meta_to_path(meta)
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir) / "temp"
        temp_path.mkdir(parents=True, exist_ok=True)
        result = task.execute(scratch_dir=temp_path, fetch=fetch, wf=wf)
        target_path.parent.mkdir(exist_ok=True, parents=True)
        result = _move_asset(result, target_path)
        logger.debug(f"Saved asset {task.asset_id} to {target_path}.")
    rss_before_kb = _current_rss_kb()
    collected = gc.collect()
    rss_after_kb = _current_rss_kb()
    logger.debug(
        f"[gc-experiment] post-task gc.collect for {task.asset_id}: "
        f"unreachable_collected={collected} rss_before_kb={rss_before_kb} rss_after_kb={rss_after_kb}"
    )
    return result


def _move_asset[A: Asset](asset: A, dst: Path) -> A:
    if isinstance(asset, (FileAsset, DirectoryAsset)):
        if dst.is_dir():
            shutil.rmtree(dst)
        return attr.evolve(
            asset,
            path=asset.path.rename(dst),
        )
    raise ValueError("Unknown Asset Type")
