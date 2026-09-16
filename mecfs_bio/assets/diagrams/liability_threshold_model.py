"""
Diagram for the Liability Threshold Model documentation: the normally
distributed latent liability across the population, cut by the threshold tau
into the unaffected majority and the shaded affected tail.
"""

from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.typst_diagram_task import TypstDiagramTask
from mecfs_bio.build_system.task.typst_source import TypstSource

_SOURCE_PATH = Path(__file__).parent / "liability_threshold_model.typ"

LIABILITY_THRESHOLD_MODEL_DIAGRAM = TypstDiagramTask(
    meta=DiagramFileMeta(AssetId("liability_threshold_model")),
    source=TypstSource(path=_SOURCE_PATH),
)
