"""File metadata for a tabular reference harmonized to a genome build."""

from pathlib import PurePath

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.read_spec import ReadSpec
from mecfs_bio.build_system.meta.reference_meta.panel_allele_frequency_columns import (
    PanelAlleleFrequencyColumns,
)


@frozen(slots=True)
class HarmonizableReferenceTableMeta(FileMeta):
    group: str
    sub_group: str
    sub_folder: PurePath
    extension: str
    id: AssetId = field(converter=AssetId)
    filename: str | None = None
    read_spec: ReadSpec | None = None
    harmonization_info: HarmonizationInfo | None = None
    # Set by allele-frequency panels: which column holds which ancestry's frequency.
    allele_frequency_columns: PanelAlleleFrequencyColumns | None = None

    def __attrs_post_init__(self):
        assert self.extension.startswith(".") or self.extension == ""

    @property
    def asset_id(self) -> AssetId:
        return self.id
