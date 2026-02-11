from typing import Mapping, Sequence

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.figures.fig_constants import FIGURE_DIRECTORY
from mecfs_bio.figures.figure_exporter import FigureExporter
from mecfs_bio.figures.figure_tasks import ALL_FIGURE_TASKS

_imo_hasher_128 = ImoHasher.with_xxhash_128()


def _runner_func(tasks: Sequence[Task]) -> Mapping[AssetId, Asset]:
    return DEFAULT_RUNNER.run(targets=list(tasks), incremental_save=True)


FIGURE_EXPORTER = FigureExporter(runner=_runner_func, tracer=_imo_hasher_128)


def generate_figures(tasks: list[Task]):
    """
    Use the build system to generate the figure assets, then copy them to the figure directory
    """
    FIGURE_EXPORTER.export(tasks, fig_dir=FIGURE_DIRECTORY)


if __name__ == "__main__":
    generate_figures(ALL_FIGURE_TASKS)
