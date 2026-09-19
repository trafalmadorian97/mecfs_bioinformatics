"""File metadata for a tabular reference oriented to a genome build (LD panel labels, the
1kg allele-frequency panel). A sibling of ReferenceFileMeta (not a subclass): generic
transformer tasks that isinstance-check ReferenceFileMeta will not capture it, so piping one
through a generic transformer fails closed. Use RenameColsTask for renames."""

from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec


@frozen
class HarmonizableReferenceTableMeta(FileMeta):
    group: str
    sub_group: str
    sub_folder: PurePath
    extension: str
    id: AssetId = field(converter=AssetId)
    filename: str | None = None
    read_spec: ReadSpec | None = None
    harmonization_info: HarmonizationInfo | None = None

    def __attrs_post_init__(self):
        assert self.extension.startswith(".") or self.extension == ""

    @property
    def asset_id(self) -> AssetId:
        return self.id
