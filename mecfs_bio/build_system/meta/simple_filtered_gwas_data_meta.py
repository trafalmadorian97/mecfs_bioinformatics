from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec


@frozen
class SimpleFilteredGWASDataMeta(FilteredGWASDataMeta):
    read_spec: ReadSpec | None = None

    @property
    def asset_id(self) -> AssetId:
        return self.id

    @classmethod
    def create(cls, asset_id: str, trait: str, project: str, sub_dir: str | PurePath):
        return cls(
            id=AssetId(asset_id),
            trait=trait,
            project=project,
            sub_dir=sub_dir,
        )
