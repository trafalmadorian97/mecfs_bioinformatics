"""
Diagram for the Linkage Disequilibrium documentation: in a non-recombining
region, two variants are correlated when they arose on the same branch of the
genealogy, regardless of the physical distance between them.

A simpler redraw of Figure 7 of the HapMap paper, dropping the SNP-position row
and the multiplicity column to keep only the same-branch / same-colour point.
"""

from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.typst_diagram_task import TypstDiagramTask
from mecfs_bio.build_system.task.typst_source import TypstSource

_SOURCE_PATH = Path(__file__).parent / "ld_same_branch_correlation.typ"

LD_SAME_BRANCH_CORRELATION_DIAGRAM = TypstDiagramTask(
    meta=DiagramFileMeta(AssetId("ld_same_branch_correlation")),
    source=TypstSource(path=_SOURCE_PATH),
)
