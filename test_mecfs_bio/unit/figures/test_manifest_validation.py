from pathlib import Path, PurePath

import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.markdown_file_meta import MarkdownFileMeta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.figures.manifest import FigureManifest
from mecfs_bio.figures.manifest_validation import (
    ManifestTaskMismatchError,
    validate_manifest_subset_of_tasks,
)


def _file_task(asset_id: str, extension: str = ".html") -> FakeTask:
    return FakeTask(
        meta=GWASPlotFileMeta(
            trait="t",
            project="p",
            extension=extension,
            id=AssetId(asset_id),
        )
    )


def _dir_task(asset_id: str) -> FakeTask:
    return FakeTask(
        meta=GWASPlotDirectoryMeta(
            trait="t",
            project="p",
            id=AssetId(asset_id),
        )
    )


def _md_task(asset_id: str) -> FakeTask:
    return FakeTask(
        meta=MarkdownFileMeta(
            id=AssetId(asset_id),
            trait="t",
            project="p",
            sub_dir=PurePath("analysis/markdown"),
        )
    )


def test_validate_accepts_exact_file_match(tmp_path: Path):
    tasks = [_file_task("plot_a"), _md_task("md_a")]
    manifest = FigureManifest(
        figures={Path("plot_a.html"): "h1", Path("md_a.mdx"): "h2"}
    )
    validate_manifest_subset_of_tasks(manifest=manifest, tasks=tasks, fig_dir=tmp_path)


def test_validate_accepts_file_under_directory_task(tmp_path: Path):
    tasks = [_dir_task("plot_dir")]
    manifest = FigureManifest(
        figures={
            Path("plot_dir/inner.html"): "h1",
            Path("plot_dir/nested/other.png"): "h2",
        }
    )
    validate_manifest_subset_of_tasks(manifest=manifest, tasks=tasks, fig_dir=tmp_path)


def test_validate_raises_on_orphan_path(tmp_path: Path):
    tasks = [_file_task("plot_a")]
    manifest = FigureManifest(
        figures={Path("plot_a.html"): "h1", Path("ghost.html"): "h2"}
    )
    with pytest.raises(ManifestTaskMismatchError) as exc_info:
        validate_manifest_subset_of_tasks(
            manifest=manifest, tasks=tasks, fig_dir=tmp_path
        )
    assert "ghost.html" in str(exc_info.value)
    assert "plot_a.html" not in str(exc_info.value)


def test_validate_does_not_treat_directory_prefix_as_file_match(tmp_path: Path):
    # "plot_dir_extra.html" must not be accepted just because a directory
    # task produces "plot_dir/" --- Path.parents only contains true
    # ancestors, so this comes for free with the Path-based check.
    tasks = [_dir_task("plot_dir")]
    manifest = FigureManifest(figures={Path("plot_dir_extra.html"): "h1"})
    with pytest.raises(ManifestTaskMismatchError):
        validate_manifest_subset_of_tasks(
            manifest=manifest, tasks=tasks, fig_dir=tmp_path
        )


def test_validate_accepts_empty_manifest(tmp_path: Path):
    validate_manifest_subset_of_tasks(
        manifest=FigureManifest.empty(), tasks=[_file_task("plot_a")], fig_dir=tmp_path
    )
