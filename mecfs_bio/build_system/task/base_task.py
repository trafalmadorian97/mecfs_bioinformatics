from abc import ABC, abstractmethod
from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.wf.base_wf import WF


class GeneratingTask(ABC):
    """
    Instructions for materializing an asset.
    """

    @property
    @abstractmethod
    def meta(self) -> Meta:
        """
        Metadata describing the target asset.
        """
        pass

    @property
    @abstractmethod
    def deps(self) -> list["Task"]:
        """
        List of tasks whose assets are needed to produce the target asset.
        """
        pass

    @property
    def asset_id(self) -> AssetId:
        return self.meta.asset_id

    @abstractmethod
    def execute(
        self,
        scratch_dir: Path,
        fetch: Fetch,
        wf: WF,
    ) -> Asset:
        """
        Materialize the target asset, using the 'fetch' callback to access asset dependencies.
        """
        pass


Task = GeneratingTask
