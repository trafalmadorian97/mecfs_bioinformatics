"""
Metadata describing tabular data from the analysis of GWAS.
"""

from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec


@frozen
class ResultTableMeta(FileMeta):
    id: str
    trait: str
    project: str
    extension: str
    read_spec: DataFrameReadSpec | None = None
    sub_dir: PurePath = PurePath("analysis")

    @property
    def asset_id(self) -> AssetId:
        return AssetId(self.id)
