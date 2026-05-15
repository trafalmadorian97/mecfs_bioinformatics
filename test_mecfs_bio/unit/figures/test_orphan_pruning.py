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


def _write_manifest(manifest_path: Path, entries: dict[str, str]) -> None:
    FigureManifest(figures=entries).save(manifest_path)


def test_find_doc_references_picks_up_md_and_mdx(tmp_path: Path):
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "a.md").write_text("![foo](../_figs/my_plot.png)\n")
    sub = docs / "sub"
    sub.mkdir()
    (sub / "b.mdx").write_text('{{ plotly_embed("docs/_figs/my_plot.png", id="x") }}\n')
    (docs / "unrelated.md").write_text("nothing here\n")

    hits = find_doc_references(rel_path="my_plot.png", docs_dir=docs)

    assert sorted(hits) == sorted([docs / "a.md", sub / "b.mdx"])


def test_find_doc_references_ignores_non_doc_files(tmp_path: Path):
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "notes.txt").write_text("references my_plot.png\n")

    assert find_doc_references(rel_path="my_plot.png", docs_dir=docs) == []


def test_prune_drops_orphan_with_no_doc_reference_and_present_file(tmp_path: Path):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    manifest_path = tmp_path / "manifest.json"

    (fig_dir / "kept.png").write_bytes(b"x")
    (fig_dir / "orphan.png").write_bytes(b"y")
    _write_manifest(manifest_path, {"kept.png": "h1", "orphan.png": "h2"})

    prune_orphan_figures(
        manifest_path=manifest_path,
        docs_dir=docs,
        fig_dir=fig_dir,
        figure_tasks=[_file_task("kept")],
    )

    assert not (fig_dir / "orphan.png").exists()
    assert (fig_dir / "kept.png").exists()
    assert FigureManifest.load(manifest_path).figures == {"kept.png": "h1"}


def test_prune_drops_orphan_when_local_file_already_missing(tmp_path: Path):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    manifest_path = tmp_path / "manifest.json"

    _write_manifest(manifest_path, {"orphan.png": "h2"})

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
    _write_manifest(manifest_path, {"orphan.png": "h2"})

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
    assert FigureManifest.load(manifest_path).figures == {"orphan.png": "h2"}


def test_prune_errors_with_missing_file_note_when_orphan_referenced_and_absent(
    tmp_path: Path,
):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "page.md").write_text("![p](../_figs/orphan.png)\n")
    manifest_path = tmp_path / "manifest.json"

    _write_manifest(manifest_path, {"orphan.png": "h2"})

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


def test_prune_is_noop_when_no_orphans(tmp_path: Path):
    fig_dir = tmp_path / "_figs"
    fig_dir.mkdir()
    docs = tmp_path / "docs"
    docs.mkdir()
    manifest_path = tmp_path / "manifest.json"

    (fig_dir / "kept.png").write_bytes(b"x")
    _write_manifest(manifest_path, {"kept.png": "h1"})

    prune_orphan_figures(
        manifest_path=manifest_path,
        docs_dir=docs,
        fig_dir=fig_dir,
        figure_tasks=[_file_task("kept")],
    )

    assert FigureManifest.load(manifest_path).figures == {"kept.png": "h1"}
    assert (fig_dir / "kept.png").exists()
