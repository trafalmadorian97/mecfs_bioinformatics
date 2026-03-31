"""
A Task that does nothing, which can be used for testing.
"""

from pathlib import Path

from attrs import Factory, frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

@frozen
class FakeTask(Task):
    meta: Meta
    deps: list[Task] = Factory(list)
    """
    For testing
    """

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        raise NotImplementedError()

    @classmethod
    def create_with_filemeta(cls, asset_id: str):
        return cls(SimpleFileMeta(AssetId(asset_id)))
