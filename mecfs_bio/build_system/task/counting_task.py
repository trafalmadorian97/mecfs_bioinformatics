"""
For testing purposes, recording the number of times that a wrapped task has been executed.
"""

from pathlib import Path

from attrs import define

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.wf.base_wf import WF


@define
class CountingTask(GeneratingTask):
    """
    For testing.  Records the number of times a task has been executed
    """

    wrapped: Task
    run_count: int = 0

    @property
    def deps(self) -> list["Task"]:
        return self.wrapped.deps

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        self.run_count += 1
        return self.wrapped.execute(scratch_dir, fetch, wf)

    @property
    def meta(self) -> Meta:
        return self.wrapped.meta
