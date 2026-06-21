"""
Metadata describing a compressed archive of analysis results.
"""

from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta


@frozen
class ResultArchiveMeta(FileMeta):
    id: AssetId
    trait: str
    project: str
    extension: str
    sub_dir: PurePath = PurePath("analysis/result_archives")

    def __attrs_post_init__(self):
        assert self.extension.startswith(".")

    @property
    def asset_id(self) -> AssetId:
        return AssetId(self.id)
