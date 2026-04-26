import hashlib
from pathlib import Path

import pytest

from mecfs_bio.figures.manifest import (
    MANIFEST_VERSION,
    FigureManifest,
    scan_figure_dir,
    sha256_of_file,
)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def test_sha256_of_file_matches_hashlib(tmp_path: Path):
    p = tmp_path / "file.bin"
    p.write_bytes(b"hello world")
    assert sha256_of_file(p) == _sha256(b"hello world")


def test_scan_figure_dir_picks_up_files_recursively(tmp_path: Path):
    (tmp_path / "a.png").write_bytes(b"AAA")
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "b.html").write_bytes(b"BBB")

    manifest = scan_figure_dir(tmp_path)
    assert manifest.figures == {
        "a.png": _sha256(b"AAA"),
        "sub/b.html": _sha256(b"BBB"),
    }


def test_manifest_load_returns_empty_when_missing(tmp_path: Path):
    assert FigureManifest.load(tmp_path / "no_such.json") == FigureManifest.empty()


def test_manifest_save_then_load_roundtrips(tmp_path: Path):
    manifest_path = tmp_path / "m.json"
    original = FigureManifest(figures={"x.png": "deadbeef", "y/z.html": "cafef00d"})
    original.save(manifest_path)
    assert FigureManifest.load(manifest_path) == original


def test_manifest_save_writes_sorted_keys(tmp_path: Path):
    manifest_path = tmp_path / "m.json"
    FigureManifest(figures={"b.png": "1", "a.png": "2"}).save(manifest_path)
    text = manifest_path.read_text()
    assert text.index('"a.png"') < text.index('"b.png"')


def test_manifest_load_rejects_unknown_version(tmp_path: Path):
    manifest_path = tmp_path / "m.json"
    manifest_path.write_text(f'{{"version": {MANIFEST_VERSION + 99}, "figures": {{}}}}')
    with pytest.raises(ValueError):
        FigureManifest.load(manifest_path)


def test_manifest_with_and_without_entry_are_pure():
    base = FigureManifest(figures={"a.png": "1"})
    added = base.with_entry("b.png", "2")
    removed = added.without_entry("a.png")
    assert base.figures == {"a.png": "1"}
    assert added.figures == {"a.png": "1", "b.png": "2"}
    assert removed.figures == {"b.png": "2"}


def test_manifest_hashes_returns_unique_set():
    manifest = FigureManifest(figures={"a.png": "h1", "b.png": "h2", "c.png": "h1"})
    assert manifest.hashes() == {"h1", "h2"}
