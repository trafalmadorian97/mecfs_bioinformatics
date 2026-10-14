"""
Diagram for the Liability Threshold Model documentation: the distribution of
liability in a case/control ascertained sample. The density above the threshold
tau is scaled up and the density below it scaled down, so that cases and
controls each make up half of the sample.
"""

from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.typst_diagram_task import TypstDiagramTask
from mecfs_bio.build_system.task.typst_source import TypstSource

_SOURCE_PATH = Path(__file__).parent / "liability_threshold_ascertained_sample.typ"

LIABILITY_THRESHOLD_ASCERTAINED_SAMPLE_DIAGRAM = TypstDiagramTask(
    meta=DiagramFileMeta(AssetId("liability_threshold_ascertained_sample_diagram")),
    source=TypstSource(path=_SOURCE_PATH),
)
