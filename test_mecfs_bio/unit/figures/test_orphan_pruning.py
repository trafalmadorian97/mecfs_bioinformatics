from pathlib import Path

import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.figures.manifest import FigureManifest
from mecfs_bio.figures.orphan_pruning import (
    OrphanReferencedInDocsError,
    find_doc_references,
    prune_orphan_figures,
)


def _file_task(asset_id: str, extension: str = ".png") -> FakeTask:
    return FakeTask(
        meta=GWASPlotFileMeta(
            trait="t",
            project="p",
            extension=extension,
            id=AssetId(asset_id),
        )
    )


def _write_manifest(manifest_path: Path, entries: dict[Path, str]) -> None:
    FigureManifest(figures=entries).save(manifest_path)


def test_find_doc_references_picks_up_md(tmp_path: Path):
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "a.md").write_text("![foo](../_figs/my_plot.png)\n")
    sub = docs / "sub"
    sub.mkdir()
    (sub / "b.md").write_text('{{ plotly_embed("docs/_figs/my_plot.png", id="x") }}\n')
    (docs / "unrelated.md").write_text("nothing here\n")

    hits = find_doc_references(rel_path=Path("my_plot.png"), docs_dir=docs)

    assert sorted(hits) == sorted([docs / "a.md", sub / "b.md"])


def test_find_doc_references_ignores_mdx_files(tmp_path: Path):
    # .mdx files in this repo are figure tables, not documentation pages,
    # so they should not gate orphan pruning.
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "table.mdx").write_text("includes my_plot.png\n")

    assert find_doc_references(rel_path=Path("my_plot.png"), docs_dir=docs) == []


def test_find_doc_references_ignores_non_doc_files(tmp_path: Path):
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "notes.txt").write_text("references my_plot.png\n")

    assert find_doc_references(rel_path=Path("my_plot.png"), docs_dir=docs) == []


def test_find_doc_references_raises_on_non_utf8_doc(tmp_path: Path):
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "bad.md").write_bytes(b"\xff\xfe not utf8")

    with pytest.raises(UnicodeDecodeError):
        find_doc_references(rel_path=Path("my_plot.png"), docs_dir=docs)


def test_prune_drops_orphan_with_no_doc_reference_and_present_file(tmp_path: Path):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    manifest_path = tmp_path / "manifest.json"

    (fig_dir / "kept.png").write_bytes(b"x")
    (fig_dir / "orphan.png").write_bytes(b"y")
    _write_manifest(manifest_path, {Path("kept.png"): "h1", Path("orphan.png"): "h2"})

    prune_orphan_figures(
        manifest_path=manifest_path,
        docs_dir=docs,
        fig_dir=fig_dir,
        figure_tasks=[_file_task("kept")],
    )

    assert not (fig_dir / "orphan.png").exists()
    assert (fig_dir / "kept.png").exists()
    assert FigureManifest.load(manifest_path).figures == {Path("kept.png"): "h1"}


def test_prune_drops_orphan_when_local_file_already_missing(tmp_path: Path):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    manifest_path = tmp_path / "manifest.json"

    _write_manifest(manifest_path, {Path("orphan.png"): "h2"})

    prune_orphan_figures(
        manifest_path=manifest_path,
        docs_dir=docs,
        fig_dir=fig_dir,
        figure_tasks=[],
    )

    assert FigureManifest.load(manifest_path).figures == {}


def test_prune_errors_when_orphan_referenced_with_present_file(tmp_path: Path):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "page.md").write_text("![p](../_figs/orphan.png)\n")
    manifest_path = tmp_path / "manifest.json"

    (fig_dir / "orphan.png").write_bytes(b"y")
    _write_manifest(manifest_path, {Path("orphan.png"): "h2"})

    with pytest.raises(OrphanReferencedInDocsError) as exc:
        prune_orphan_figures(
            manifest_path=manifest_path,
            docs_dir=docs,
            fig_dir=fig_dir,
            figure_tasks=[],
        )

    msg = str(exc.value)
    assert "orphan.png" in msg
    assert "page.md" in msg
    assert "local file is also missing" not in msg
    # Disk untouched on error
    assert (fig_dir / "orphan.png").exists()
    assert FigureManifest.load(manifest_path).figures == {Path("orphan.png"): "h2"}


def test_prune_errors_with_missing_file_note_when_orphan_referenced_and_absent(
    tmp_path: Path,
):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "page.md").write_text("![p](../_figs/orphan.png)\n")
    manifest_path = tmp_path / "manifest.json"

    _write_manifest(manifest_path, {Path("orphan.png"): "h2"})

    with pytest.raises(OrphanReferencedInDocsError) as exc:
        prune_orphan_figures(
            manifest_path=manifest_path,
            docs_dir=docs,
            fig_dir=fig_dir,
            figure_tasks=[],
        )

    msg = str(exc.value)
    assert "orphan.png" in msg
    assert "page.md" in msg
    assert "local file is also missing" in msg


def test_prune_detects_orphan_file_not_yet_in_manifest(tmp_path: Path):
    # A leftover file under fig_dir whose task was removed from
    # ALL_FIGURE_TASKS may never have been written to the manifest. The
    # manifest-merge inside push_figures would still pick it up via its
    # scan, so pruning has to catch it from the on-disk side too.
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    manifest_path = tmp_path / "manifest.json"

    (fig_dir / "leftover.png").write_bytes(b"x")
    _write_manifest(manifest_path, {})

    prune_orphan_figures(
        manifest_path=manifest_path,
        docs_dir=docs,
        fig_dir=fig_dir,
        figure_tasks=[],
    )

    assert not (fig_dir / "leftover.png").exists()
    assert FigureManifest.load(manifest_path).figures == {}


def test_prune_errors_on_referenced_leftover_file_not_in_manifest(tmp_path: Path):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "page.md").write_text("![p](../_figs/leftover.png)\n")
    manifest_path = tmp_path / "manifest.json"

    (fig_dir / "leftover.png").write_bytes(b"x")
    _write_manifest(manifest_path, {})

    with pytest.raises(OrphanReferencedInDocsError) as exc:
        prune_orphan_figures(
            manifest_path=manifest_path,
            docs_dir=docs,
            fig_dir=fig_dir,
            figure_tasks=[],
        )
    assert "leftover.png" in str(exc.value)
    assert (fig_dir / "leftover.png").exists()


def test_prune_is_noop_when_no_orphans(tmp_path: Path):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    manifest_path = tmp_path / "manifest.json"

    (fig_dir / "kept.png").write_bytes(b"x")
    _write_manifest(manifest_path, {Path("kept.png"): "h1"})

    prune_orphan_figures(
        manifest_path=manifest_path,
        docs_dir=docs,
        fig_dir=fig_dir,
        figure_tasks=[_file_task("kept")],
    )

    assert FigureManifest.load(manifest_path).figures == {Path("kept.png"): "h1"}
    assert (fig_dir / "kept.png").exists()
