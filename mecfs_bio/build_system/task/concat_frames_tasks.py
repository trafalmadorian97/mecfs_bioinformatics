from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class ConcatFramesTask(Task):
    @property
    def meta(self) -> Meta:
        pass

    @property
    def deps(self) -> list["Task"]:
        pass

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        pass