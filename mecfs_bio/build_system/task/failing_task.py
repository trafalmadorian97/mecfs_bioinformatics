"""
A task that always fails, which can be used for testing.
"""

from pathlib import Path
from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

@frozen
class FailingTask(Task):
    """
    A task that always fails.
    Intended to be used for testing.`
    """

    meta: Meta
    deps: list[Task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        raise ValueError("Task intentionally failed.")
