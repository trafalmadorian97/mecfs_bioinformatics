"""
Task that renders a matplotlib figure into an SVG diagram asset.
"""

from pathlib import Path
from typing import Protocol

import matplotlib
from attrs import frozen
from matplotlib.figure import Figure

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.diagram_file_meta import (
    DIAGRAM_EXTENSION,
    DiagramFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.wf.base_wf import WF


class FigureBuilder(Protocol):
    """
    Anything that can construct a complete matplotlib figure from scratch.

    build should create the figure with matplotlib.figure.Figure rather than
    pyplot, so that no figure is left registered in pyplot's global state.
    """

    def build(self) -> Figure: ...


@frozen(slots=True)
class MatplotlibDiagramTask(GeneratingTask):
    """
    Builds a matplotlib figure and saves it as an SVG diagram.

    Text is written as paths (svg.fonttype = path), so the SVG renders
    identically without the fonts installed.

    By default, matplotlib's SVG output differs from run to run even when the
    figure is unchanged, for two reasons. First, it records the creation date
    in the file's metadata. Second, the ids it gives to reusable elements
    (glyphs, clip paths, markers) are hashes that mix in a random value drawn
    afresh on every save. Byte-level changes like these would make an unchanged
    figure look new to anything that compares file contents, such as the
    figures manifest. This task omits the date (metadata Date = None) and
    replaces the random value with the asset id (svg.hashsalt), so rebuilding
    unchanged code yields byte-identical output. Using the asset id rather than
    one shared constant also keeps ids distinct between different figures, in
    case several SVGs are inlined into one HTML page, where ids must be unique.
    """

    meta: DiagramFileMeta
    figure: FigureBuilder

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        out_path = scratch_dir / (self.meta.asset_id + DIAGRAM_EXTENSION)
        fig = self.figure.build()
        with matplotlib.rc_context(
            {"svg.fonttype": "path", "svg.hashsalt": self.meta.asset_id}
        ):
            fig.savefig(out_path, format="svg", metadata={"Date": None})
        return FileAsset(out_path)
