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
    extension: str
    id: AssetId = field(converter=AssetId)
    filename: str | None = None
    read_spec: ReadSpec | None = None

    def __attrs_post_init__(self):
        assert self.extension.startswith(".")

    @property
    def asset_id(self) -> AssetId:
        return self.id
