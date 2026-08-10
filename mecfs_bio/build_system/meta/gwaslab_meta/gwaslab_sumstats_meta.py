from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta


@frozen
class GWASLabSumStatsMeta(FileMeta):
    """
    Metadata describing a pickled gwaslab summary stats object
    """

    @property
    def asset_id(self) -> AssetId:
        return self.id

    id: AssetId
    trait: str
    project: str
    sub_dir: str = "gwaslab_sumstats"
