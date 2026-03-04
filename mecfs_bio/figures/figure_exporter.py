import shutil
from pathlib import Path
from typing import Callable, Mapping, Sequence

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.base_tracer import (
    Tracer,
)
from mecfs_bio.build_system.task.base_task import Task

logger = structlog.get_logger()


@frozen
class FigureExporter:
    """
    Responsible for invoking the build system to generate Assets corresponding to figures, then copying those assets to the figure directory.
    """

    runner: Callable[[Sequence[Task]], Mapping[AssetId, Asset]]
    tracer: Tracer

    def export(self, to_export: Sequence[Task], fig_dir: Path) -> None:
        fig_dir.mkdir(parents=True, exist_ok=True)
        meta_list = []
        for task in to_export:
            assert isinstance(task.meta, (GWASPlotFileMeta, GWASPlotDirectoryMeta))
            meta_list.append(task.meta)
        result = self.runner(
            to_export,
        )
        for task, meta in zip(to_export, meta_list):
            asset = result[task.asset_id]
            if isinstance(asset, FileAsset):
                src = asset.path
                assert isinstance(meta, GWASPlotFileMeta)
                dst = get_fig_file_path(meta=meta, fig_dir=fig_dir)
                shutil.copy(src, dst)
                logger.debug(f"Figure asset {task.asset_id} copied to {dst}.")
            elif isinstance(asset, DirectoryAsset):
                src = asset.path
                assert isinstance(meta, GWASPlotDirectoryMeta)
                dst = get_fig_dir_meta(meta=meta, fig_dir=fig_dir)
                shutil.copytree(src, dst)
                logger.debug(f"Directory figure asset {task.asset_id} copied to {dst}.")


def get_fig_file_path(meta: GWASPlotFileMeta, fig_dir: Path) -> Path:
    return fig_dir / (meta.asset_id + meta.extension)


def get_fig_dir_meta(meta: GWASPlotDirectoryMeta, fig_dir: Path) -> Path:
    return fig_dir / meta.id
