"""
Diagram for the colocalization (coloc.abf) documentation: the configuration
space of causal-variant pairs (S_1, S_2), each cell coloured by the hypothesis
it belongs to.

S_k in {0, 1, ..., Q} is the causal variant of trait k in a locus of Q variants,
with S_k = 0 meaning trait k has none. The (Q + 1)^2 configurations partition
into H_0 (neither trait), H_1 (trait 1 only), H_2 (trait 2 only), H_3 (both
traits, distinct variants) and H_4 (both traits, one shared variant, the
diagonal).
"""

from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.typst_diagram_task import TypstDiagramTask
from mecfs_bio.build_system.task.typst_source import TypstSource

_SOURCE_PATH = Path(__file__).parent / "coloc_configuration_grid.typ"

COLOC_CONFIGURATION_GRID_DIAGRAM = TypstDiagramTask(
    meta=DiagramFileMeta(AssetId("coloc_configuration_grid")),
    source=TypstSource(path=_SOURCE_PATH),
)
