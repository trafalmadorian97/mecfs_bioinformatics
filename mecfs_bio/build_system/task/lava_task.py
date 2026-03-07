from pathlib import Path
from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import PhenotypeInfo
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF



@frozen
class LavaDataSource:
    """
    A source Task providing tabular Gwas summary statistics data
    """

    task: Task
    alias: str
    pipe: DataProcessingPipe = IdentityPipe()
    sample_info: PhenotypeInfo | None = None

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id

@frozen
class LavaTask(Task):
    _meta: Meta
    sources: Sequence[LavaDataSource]
    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [item.task for item in self.sources]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        pass