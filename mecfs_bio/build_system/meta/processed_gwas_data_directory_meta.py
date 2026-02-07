"""
Metadata describing a directory containing plots from the analysis of GWAS data.
"""

from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta


@frozen
class ProcessedGwasDataDirectoryMeta(DirMeta):
    @property
    def asset_id(self) -> AssetId:
        return self.id

    id: AssetId
    trait: str
    project: str
    sub_dir: PurePath
