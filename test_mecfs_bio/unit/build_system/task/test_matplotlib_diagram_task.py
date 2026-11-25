from pathlib import Path

from attrs import frozen
from matplotlib.figure import Figure

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.matplotlib_diagram_task import MatplotlibDiagramTask
from mecfs_bio.build_system.wf.base_wf import make_wf


@frozen(slots=True)
class _LineFigure:
    def build(self) -> Figure:
        fig = Figure(figsize=(1.0, 1.0))
        ax = fig.add_subplot()
        ax.plot([0.0, 1.0], [0.0, 1.0])
        ax.set_title("smoke")
        return fig


def _no_fetch(asset_id: AssetId) -> Asset:
    raise ValueError("MatplotlibDiagramTask has no dependencies")


def _render(scratch: Path) -> FileAsset:
    scratch.mkdir()
    task = MatplotlibDiagramTask(
        meta=DiagramFileMeta(AssetId("smoke")), figure=_LineFigure()
    )
    result = task.execute(scratch_dir=scratch, fetch=_no_fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return result


def test_renders_svg(tmp_path: Path) -> None:
    result = _render(tmp_path / "scratch")
    head = result.path.read_bytes()[:256].lstrip()
    assert head.startswith(b"<svg") or head.startswith(b"<?xml")


def test_rerendering_is_byte_identical(tmp_path: Path) -> None:
    first = _render(tmp_path / "first")
    second = _render(tmp_path / "second")
    assert first.path.read_bytes() == second.path.read_bytes()
