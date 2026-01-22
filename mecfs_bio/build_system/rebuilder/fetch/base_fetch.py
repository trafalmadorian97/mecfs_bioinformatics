from abc import abstractmethod
from typing import Protocol

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId


class Fetch(Protocol):
    """
    An interface for materializing or retrieving assets, given their ids.
    """

    @abstractmethod
    def __call__(self, asset_id: AssetId) -> Asset:
        pass
