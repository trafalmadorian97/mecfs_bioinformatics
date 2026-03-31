"""
Metadata describing GWAS summary statistics.
"""

from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec


@frozen
class GWASSummaryDataFileMeta(FileMeta):
    @property
    def asset_id(self) -> AssetId:
        return self.id

    id: AssetId = field(converter=AssetId)
    trait: str
    project: str
    sub_dir: str
    project_path: PurePath | None
    read_spec: ReadSpec | None = None
    extension: str | None = None

