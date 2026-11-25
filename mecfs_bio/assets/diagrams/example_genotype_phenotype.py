"""
Example diagram: a genotype node pointing to a phenotype node. Establishes the
diagram-source convention; safe to replace with a real documentation diagram.
"""

from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.typst_diagram_task import TypstDiagramTask
from mecfs_bio.build_system.task.typst_source import TypstSource

_SOURCE_PATH = Path(__file__).parent / "example_genotype_phenotype.typ"

EXAMPLE_GENOTYPE_PHENOTYPE_DIAGRAM = TypstDiagramTask(
    meta=DiagramFileMeta(AssetId("example_genotype_phenotype")),
    source=TypstSource(path=_SOURCE_PATH),
)
