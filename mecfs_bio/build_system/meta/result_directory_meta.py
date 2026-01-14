"""
Metadata describing a directory containing multiple analysis result files

"""
from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta


@frozen
class ResultDirectoryMeta(DirMeta):
    _asset_id: str
    trait: str
    project: str
    sub_dir: PurePath = PurePath("analysis")

    @property
    def asset_id(self) -> AssetId:
        return AssetId(self._asset_id)

