from pathlib import Path, PurePath
from typing import Sequence

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_manhattan_plot_meta import (
    GWASLabManhattanQQPlotMeta,
)
from mecfs_bio.build_system.meta.markdown_file_meta import MarkdownFileMeta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.figures.figure_exporter import AbstractFigureExporter
from mecfs_bio.figures.key_scripts.generate_new_figures import generate_new_figures


class _RecordingExporter(AbstractFigureExporter):
    """Records the tasks passed to export() instead of doing real work."""

    def __init__(self) -> None:
        self.calls: list[tuple[tuple[Task, ...], Path]] = []

    def export(self, to_export: Sequence[Task], fig_dir: Path) -> None:
        self.calls.append((tuple(to_export), fig_dir))


def test_generate_new_figures_only_passes_tasks_with_missing_outputs(
    tmp_path: Path,
):
    fig_dir = tmp_path

    existing_plot_task = FakeTask(
        meta=GWASPlotFileMeta(
            trait="trait_a",
            project="proj_a",
            extension=".html",
            id=AssetId("existing_plot"),
        )
    )
    missing_plot_task = FakeTask(
        meta=GWASPlotFileMeta(
            trait="trait_b",
            project="proj_b",
            extension=".html",
            id=AssetId("missing_plot"),
        )
    )
    existing_dir_task = FakeTask(
        meta=GWASPlotDirectoryMeta(
            trait="trait_c",
            project="proj_c",
            id=AssetId("existing_dir"),
        )
    )
    missing_dir_task = FakeTask(
        meta=GWASPlotDirectoryMeta(
            trait="trait_d",
            project="proj_d",
            id=AssetId("missing_dir"),
        )
    )
    existing_md_task = FakeTask(
        meta=MarkdownFileMeta(
            id=AssetId("existing_md"),
            trait="trait_e",
            project="proj_e",
            sub_dir=PurePath("analysis/markdown"),
        )
    )
    missing_md_task = FakeTask(
        meta=MarkdownFileMeta(
            id=AssetId("missing_md"),
            trait="trait_f",
            project="proj_f",
            sub_dir=PurePath("analysis/markdown"),
        )
    )
    existing_manhattan_task = FakeTask(
        meta=GWASLabManhattanQQPlotMeta(
            trait="trait_g",
            project="proj_g",
            id=AssetId("existing_manhattan"),
        )
    )
    missing_manhattan_task = FakeTask(
        meta=GWASLabManhattanQQPlotMeta(
            trait="trait_h",
            project="proj_h",
            id=AssetId("missing_manhattan"),
        )
    )

    (fig_dir / "existing_plot.html").write_text("placeholder")
    (fig_dir / "existing_dir").mkdir()
    (fig_dir / "existing_md.mdx").write_text("placeholder")
    (fig_dir / "existing_manhattan.png").write_text("placeholder")

    all_tasks: list[Task] = [
        existing_plot_task,
        missing_plot_task,
        existing_dir_task,
        missing_dir_task,
        existing_md_task,
        missing_md_task,
        existing_manhattan_task,
        missing_manhattan_task,
    ]

    exporter = _RecordingExporter()
    generate_new_figures(all_figure_tasks=all_tasks, fig_dir=fig_dir, exporter=exporter)

    assert len(exporter.calls) == 1
    exported_tasks, called_fig_dir = exporter.calls[0]
    assert called_fig_dir == fig_dir
    assert {task.asset_id for task in exported_tasks} == {
        missing_plot_task.asset_id,
        missing_dir_task.asset_id,
        missing_md_task.asset_id,
        missing_manhattan_task.asset_id,
    }
