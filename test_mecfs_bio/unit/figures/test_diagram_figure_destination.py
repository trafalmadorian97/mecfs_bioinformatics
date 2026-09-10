from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import (
    DIAGRAM_EXTENSION,
    DiagramFileMeta,
)
from mecfs_bio.figures.figure_exporter import get_figure_destination


def test_diagram_lands_flat_in_fig_dir(tmp_path: Path):
    meta = DiagramFileMeta(AssetId("concept"))
    assert get_figure_destination(meta=meta, fig_dir=tmp_path) == tmp_path / (
        "concept" + DIAGRAM_EXTENSION
    )
