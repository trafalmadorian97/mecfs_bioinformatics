from pathlib import Path

import pytest

from mecfs_bio.build_system.task.typst_source import TypstSource


def test_rejects_both_path_and_code(tmp_path: Path):
    f = tmp_path / "d.typ"
    f.write_text("x")
    with pytest.raises(AssertionError):
        TypstSource(path=f, code="x")


def test_rejects_neither():
    with pytest.raises(AssertionError):
        TypstSource()


def test_rejects_missing_file(tmp_path: Path):
    with pytest.raises(AssertionError):
        TypstSource(path=tmp_path / "does_not_exist.typ")


def test_resolves_code_to_scratch_file(tmp_path: Path):
    src = TypstSource(code="hello typst")
    resolved = src.resolve_to_path(tmp_path)
    assert resolved.read_text() == "hello typst"


def test_returns_existing_path_unchanged(tmp_path: Path):
    f = tmp_path / "d.typ"
    f.write_text("content")
    src = TypstSource(path=f)
    assert src.resolve_to_path(tmp_path) == f
