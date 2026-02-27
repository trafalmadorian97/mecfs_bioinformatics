"""
Base classes for the two types of metadata: FileMeta and DirMeta.
"""

from abc import ABC, abstractmethod

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec


class FileMeta(ABC):
    """
    Metadata describing a file asset that either currently exists or does not exist but can be materialized.
    """

    @property
    @abstractmethod
    def asset_id(self) -> AssetId:
        """
        A uniquely identifying ID for the asset
        """
        pass

    def read_spec(self) -> ReadSpec | None:
        """
        A specifier describing how the data should be read
        """
        return None


class DirMeta(ABC):
    """
    Metadata describing a directory asset that either currently exists or does not exist but can be materialized.
    """

    @property
    @abstractmethod
    def asset_id(self) -> AssetId:
        pass
