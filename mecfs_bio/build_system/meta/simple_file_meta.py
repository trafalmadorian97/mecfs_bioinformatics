"""
Metadata describing a miscellaneous file.
"""

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec


@frozen
class SimpleFileMeta(FileMeta):
    short_id: AssetId = field(converter=AssetId)
    _read_spec: ReadSpec | None = None

    @property
    def asset_id(self) -> AssetId:
        return self.short_id

    def read_spec(self) -> ReadSpec | None:
        return self._read_spec

    @classmethod
    def create(cls, asset_id: str):
        return cls(
            short_id=AssetId(asset_id),
        )
