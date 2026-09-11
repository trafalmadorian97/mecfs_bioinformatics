"""
Metadata describing a single SVG diagram asset generated from Typst source.
"""

from attrs import field, frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import FileMeta

DIAGRAM_EXTENSION = ".svg"


@frozen
class DiagramFileMeta(FileMeta):
    """
    Metadata describing a single SVG diagram.

    A diagram is a conceptual illustration for the documentation, not tied to
    any GWAS, so unlike GWASPlotFileMeta it carries no trait or project. The
    file extension is always .svg.
    """

    id: AssetId = field(converter=AssetId)

    @property
    def asset_id(self) -> AssetId:
        return self.id
