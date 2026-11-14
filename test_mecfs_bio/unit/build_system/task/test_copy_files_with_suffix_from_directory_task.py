from pathlib import Path, PurePath

import pytest

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.task.copy_files_with_suffix_from_directory_task import (
    CopyFilesWithSuffixFromDirectoryTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf


def _source_task() -> FakeTask:
    return FakeTask(
        ResultDirectoryMeta(
            id=AssetId("plots"),
            trait="synthetic",
            project="test",
            sub_dir=PurePath("analysis"),
        )
    )


def _copy(tmp_path: Path, filenames: list[str], suffix: str) -> Path:
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    for name in filenames:
        (source_dir / name).write_text(name)
    task = CopyFilesWithSuffixFromDirectoryTask.create(
        asset_id="copied", source_directory_task=_source_task(), suffix=suffix
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch,
        fetch=lambda asset_id: {AssetId("plots"): DirectoryAsset(source_dir)}[asset_id],
        wf=make_wf(),
    )
    assert isinstance(result, DirectoryAsset)
    return result.path


def test_copies_only_files_with_suffix(tmp_path: Path) -> None:
    out = _copy(
        tmp_path, ["a.svg", "b.svg", "a.png", "index.parquet", "c.svgz"], ".svg"
    )
    assert sorted(p.name for p in out.iterdir()) == ["a.svg", "b.svg"]
    assert (out / "a.svg").read_text() == "a.svg"


def test_no_matching_files_gives_empty_directory(tmp_path: Path) -> None:
    out = _copy(tmp_path, ["index.parquet"], ".png")
    assert list(out.iterdir()) == []


def test_suffix_must_start_with_dot() -> None:
    with pytest.raises(AssertionError):
        CopyFilesWithSuffixFromDirectoryTask.create(
            asset_id="copied", source_directory_task=_source_task(), suffix="svg"
        )
