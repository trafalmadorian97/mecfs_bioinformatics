import shutil
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Callable, Mapping, Sequence

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_manhattan_plot_meta import (
    GWASLabManhattanQQPlotMeta,
)
from mecfs_bio.build_system.meta.markdown_file_meta import MarkdownFileMeta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.base_tracer import (
    Tracer,
)
from mecfs_bio.build_system.task.base_task import Task

logger = structlog.get_logger()

ValidFigureMeta = (
    GWASPlotFileMeta
    | GWASPlotDirectoryMeta
    | MarkdownFileMeta
    | GWASLabManhattanQQPlotMeta
)


class AbstractFigureExporter(ABC):
    @abstractmethod
    def export(self, to_export: Sequence[Task], fig_dir: Path) -> None:
        pass


@frozen
class FigureExporter(AbstractFigureExporter):
    """
    Responsible for invoking the build system to generate Assets corresponding to figures, then copying those assets to the figure directory.
    """

    runner: Callable[[Sequence[Task]], Mapping[AssetId, Asset]]
    tracer: Tracer

    def export(self, to_export: Sequence[Task], fig_dir: Path) -> None:
        fig_dir.mkdir(parents=True, exist_ok=True)
        meta_list = []
        for task in to_export:
            assert isinstance(task.meta, ValidFigureMeta)
            meta_list.append(task.meta)
        result = self.runner(
            to_export,
        )
        for task, meta in zip(to_export, meta_list):
            asset = result[task.asset_id]
            dst = get_figure_destination(meta=meta, fig_dir=fig_dir)
            if isinstance(meta, GWASPlotDirectoryMeta):
                assert isinstance(asset, DirectoryAsset)
                shutil.copytree(asset.path, dst, dirs_exist_ok=True)
                logger.debug(f"Directory figure asset {task.asset_id} copied to {dst}.")
            else:
                assert isinstance(asset, FileAsset)
                shutil.copy(asset.path, dst)
                logger.debug(f"Figure asset {task.asset_id} copied to {dst}.")


def get_fig_file_path(
    meta: GWASPlotFileMeta | GWASLabManhattanQQPlotMeta, fig_dir: Path
) -> Path:
    if isinstance(meta, GWASLabManhattanQQPlotMeta):
        return fig_dir / (str(meta.asset_id) + ".png")
    return fig_dir / (meta.asset_id + meta.extension)


def get_md_fig_file_path(meta: MarkdownFileMeta, fig_dir: Path) -> Path:
    return (
        fig_dir / (meta.asset_id + ".mdx")
    )  # Use .mdx instead of .md to prevent mkdocs from auto-generating a documentation page


def get_fig_dir_meta(meta: GWASPlotDirectoryMeta, fig_dir: Path) -> Path:
    return fig_dir / meta.id


def get_figure_destination(meta: ValidFigureMeta, fig_dir: Path) -> Path:
    """
    Single source of truth for "where in the figure directory does the figure
    produced by a task with this meta land?".

    For file-emitting metas this is an exact file path; for
    GWASPlotDirectoryMeta this is the destination directory (the task's
    output files all live underneath it).
    """
    if isinstance(meta, GWASPlotFileMeta | GWASLabManhattanQQPlotMeta):
        return get_fig_file_path(meta=meta, fig_dir=fig_dir)
    if isinstance(meta, GWASPlotDirectoryMeta):
        return get_fig_dir_meta(meta=meta, fig_dir=fig_dir)
    if isinstance(meta, MarkdownFileMeta):
        return get_md_fig_file_path(meta=meta, fig_dir=fig_dir)
    raise ValueError(f"Unknown meta type {type(meta)}")
