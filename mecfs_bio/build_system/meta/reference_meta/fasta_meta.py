"""Directory metadata for an indexed genome FASTA, tagged with its genome build."""

from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild


@frozen(slots=True)
class FASTAMeta(DirMeta):
    group: str
    sub_group: str
    sub_folder: PurePath
    build: GenomeBuild
    id: AssetId = field(converter=AssetId)
    dirname: str | None = None

    @property
    def asset_id(self) -> AssetId:
        return self.id
