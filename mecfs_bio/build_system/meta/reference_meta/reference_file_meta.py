from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec


@frozen
class ReferenceFileMeta(FileMeta):
    group: str
    sub_group: str
    sub_folder: PurePath
    _asset_id: AssetId = field(converter=AssetId)
    filename: str | None = None
    extension: str = ".zip"
    _read_spec: ReadSpec | None = None

    def __attrs_post_init__(self):
        assert self.extension.startswith(".")

    def read_spec(self) -> ReadSpec | None:
        return self._read_spec

    @property
    def asset_id(self) -> AssetId:
        return self._asset_id
