import shutil
import tempfile
from pathlib import Path
from typing import Callable

import attr
import structlog
from loguru import logger

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


def noop_post_execute_hook(task: Task) -> None:
    """Default post-execute hook: do nothing."""


def sandboxed_execute(
    task: Task,
    meta_to_path: MetaToPath,
    wf: WF,
    fetch: Fetch,
    post_execute: Callable[[Task], None] = noop_post_execute_hook,
) -> Asset:
    """
    Execute a task in a temporary directory, then move its result to the final location
    determined by the task's metadata.

    `post_execute` is invoked with the task after the temp directory is cleaned up and
    the asset has been moved into place. The default is a no-op.
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
    post_execute(task)
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
