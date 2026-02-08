from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta


@frozen
class GWASPlotFileMeta(FileMeta):
    trait: str
    project: str
    extension: str
    id: AssetId = field(converter=AssetId)
    sub_dir: PurePath = PurePath("analysis/plots")

    def __attrs_post_init__(self):
        assert self.extension.startswith(".")

    @property
    def asset_id(self) -> AssetId:
        return self.id
