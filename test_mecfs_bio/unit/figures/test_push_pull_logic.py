"""
Unit tests for the pure logic of push_figures and pull_figures (no network).
"""

from pathlib import Path

from mecfs_bio.figures.key_scripts.pull_figures import _prune_unmanifested
from mecfs_bio.figures.key_scripts.push_figures import _merge_manifests
from mecfs_bio.figures.manifest import FigureManifest


def test_merge_manifests_local_wins_for_overlapping_paths():
    old = FigureManifest(figures={Path("a.png"): "h_old"})
    local = FigureManifest(figures={Path("a.png"): "h_new"})
    merged = _merge_manifests(old=old, local=local, prune=False)
    assert merged.figures == {Path("a.png"): "h_new"}


def test_merge_manifests_keeps_old_only_paths_when_not_pruning():
    old = FigureManifest(figures={Path("only_in_old.png"): "h_old"})
    local = FigureManifest(figures={Path("only_in_local.png"): "h_local"})
    merged = _merge_manifests(old=old, local=local, prune=False)
    assert merged.figures == {
        Path("only_in_old.png"): "h_old",
        Path("only_in_local.png"): "h_local",
    }


def test_merge_manifests_drops_old_only_paths_when_pruning():
    old = FigureManifest(figures={Path("only_in_old.png"): "h_old"})
    local = FigureManifest(figures={Path("only_in_local.png"): "h_local"})
    merged = _merge_manifests(old=old, local=local, prune=True)
    assert merged.figures == {Path("only_in_local.png"): "h_local"}


def test_prune_unmanifested_deletes_files_not_in_manifest(tmp_path: Path):
    (tmp_path / "kept.png").write_bytes(b"x")
    (tmp_path / "stale.png").write_bytes(b"y")
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "also_stale.html").write_bytes(b"z")

    manifest = FigureManifest(figures={Path("kept.png"): "irrelevant"})
    _prune_unmanifested(fig_dir=tmp_path, manifest=manifest)

    assert (tmp_path / "kept.png").exists()
    assert not (tmp_path / "stale.png").exists()
    assert not (sub / "also_stale.html").exists()
