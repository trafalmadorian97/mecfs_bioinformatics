from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta

# from mecfs_bio.build_system.meta.file_meta import FileMeta


@frozen
class SimpleDirectoryMeta(DirMeta):
    directory_short_id: AssetId = field(converter=AssetId)

    @property
    def asset_id(self) -> AssetId:
        return self.directory_short_id
