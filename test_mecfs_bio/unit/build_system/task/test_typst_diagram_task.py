from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.diagram_file_meta import DiagramFileMeta
from mecfs_bio.build_system.task.typst_diagram_task import TypstDiagramTask
from mecfs_bio.build_system.task.typst_source import TypstSource
from mecfs_bio.build_system.wf.base_wf import make_wf

# A single-page Typst document with visible content, no @preview imports, so the
# render is hermetic (no package-registry fetch) and fast.
_SMOKE_TYPST = (
    "#set page(width: auto, height: auto, margin: 4pt)\n"
    "#circle(radius: 8pt)\n"
)


def _no_fetch(asset_id: AssetId) -> Asset:
    raise ValueError("TypstDiagramTask has no dependencies")


def test_renders_svg_from_inline_code(tmp_path: Path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    task = TypstDiagramTask(
        meta=DiagramFileMeta(AssetId("smoke")),
        source=TypstSource(code=_SMOKE_TYPST),
    )
    result = task.execute(scratch_dir=scratch, fetch=_no_fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    head = result.path.read_bytes()[:256].lstrip()
    assert head.startswith(b"<svg") or head.startswith(b"<?xml")


def test_renders_svg_from_file_source(tmp_path: Path):
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    typ_file = tmp_path / "diagram.typ"
    typ_file.write_text(_SMOKE_TYPST)
    task = TypstDiagramTask(
        meta=DiagramFileMeta(AssetId("smoke_file")),
        source=TypstSource(path=typ_file),
    )
    result = task.execute(scratch_dir=scratch, fetch=_no_fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    head = result.path.read_bytes()[:256].lstrip()
    assert head.startswith(b"<svg") or head.startswith(b"<?xml")
