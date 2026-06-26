import errno
import os
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
    if not isinstance(asset, (FileAsset, DirectoryAsset)):
        raise ValueError("Unknown Asset Type")
    if dst.is_dir():
        shutil.rmtree(dst)
    try:
        new_path = asset.path.rename(dst)
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise
        new_path = _cross_device_move(asset.path, dst)
    return attr.evolve(asset, path=new_path)


def _cross_device_move(src: Path, dst: Path) -> Path:
    """
    Move src to dst when they live on different filesystems, where os.rename fails with
    EXDEV.  The asset is first copied to a temporary sibling on dst's filesystem and then
    atomically swapped into place with os.replace, so the canonical asset path never holds
    a partially written asset.  The staging copy is removed if the copy itself fails.

    src is left in place; it lives in the caller's scratch directory, which is cleaned up
    separately.
    """
    staging = dst.parent / f".{dst.name}.partial-{os.getpid()}"
    try:
        _copy_without_metadata(src, staging)
    except BaseException:
        if staging.is_dir():
            shutil.rmtree(staging, ignore_errors=True)
        else:
            staging.unlink(missing_ok=True)
        raise
    os.replace(staging, dst)
    return dst


def _copy_without_metadata(src: Path, dst: Path) -> None:
    """
    Copy a file or directory tree's contents without replicating source timestamps or
    permissions.  Asset integrity is tracked by the build system's trace, not by
    filesystem metadata, and some asset stores live on mounts (for example Windows drives
    mounted in WSL via DrvFs) that reject utime/chmod, which makes shutil.copy2 and
    shutil.copytree fail in copystat after the data has already been copied.  Copying data
    only keeps cross-device moves working on those filesystems.
    """
    if src.is_dir():
        for dirpath, _dirnames, filenames in os.walk(src):
            rel = Path(dirpath).relative_to(src)
            (dst / rel).mkdir(parents=True, exist_ok=True)
            for name in filenames:
                shutil.copyfile(Path(dirpath) / name, dst / rel / name)
    else:
        shutil.copyfile(src, dst)
