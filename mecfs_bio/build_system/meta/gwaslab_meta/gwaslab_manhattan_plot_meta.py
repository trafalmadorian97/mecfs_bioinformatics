from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta


@frozen
class GWASLabManhattanQQPlotMeta(FileMeta):
    trait: str
    project: str
    id: AssetId = field(converter=AssetId)
    sub_dir: PurePath = PurePath("analysis/manhattan_plot")

    @property
    def asset_id(self) -> AssetId:
        return self.id
