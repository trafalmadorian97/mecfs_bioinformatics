"""
Metadata describing a markdown file
"""

from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec


@frozen
class MarkdownFileMeta(FileMeta):
    @property
    def asset_id(self) -> AssetId:
        return self.id

    id: AssetId
    trait: str
    project: str
    sub_dir: str | PurePath

    def read_spec(self) -> ReadSpec | None:
        return None
