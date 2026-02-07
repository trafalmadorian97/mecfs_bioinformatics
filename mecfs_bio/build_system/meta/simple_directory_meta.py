"""
Metadata describing a miscellaneous directory.
"""

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta


@frozen
class SimpleDirectoryMeta(DirMeta):
    id: AssetId = field(converter=AssetId)

    @property
    def asset_id(self) -> AssetId:
        return self.id
