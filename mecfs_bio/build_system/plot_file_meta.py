from pathlib import PurePath

from attrs import frozen, field

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta


@frozen
class GWASPlotFileMeta(FileMeta):
    trait: str
    project: str
    id: AssetId = field(converter=AssetId)
    sub_dir: PurePath = PurePath("analysis/plots")

    @property
    def asset_id(self) -> AssetId:
        return self.id
