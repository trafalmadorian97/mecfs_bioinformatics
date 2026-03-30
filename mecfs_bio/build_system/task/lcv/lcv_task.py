"""
Task to apply the Latent Causal Variable technique of O'Connor and Price to attempt
to estimate the causal director between two genetically correlated traits.
The key output value is GCP: Genetic Causality Proportion.

Citation:
O’Connor, Luke J., and Alkes L. Price. "Distinguishing genetic correlation
from causation across 52 diseases and complex traits." Nature genetics 50.12 (2018): 1728-1734.

"""
from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class LCVTask(Task):
    """
    Task to apply the Latent Causal Variable technique of O'Connor and Price to attempt
    to estimate the causal director between two genetically correlated traits.
    The key output value is GCP: Genetic Causality Proportion.

    """
    _meta: Meta
    trait_1_data: Task
    trait_2_data: Task

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.trait_1_data, self.trait_2_data]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        trait_1_asset = fetch(self.trait_1_data.asset_id)
        trait_2_asset = fetch(self.trait_2_data.asset_id)
        #wip