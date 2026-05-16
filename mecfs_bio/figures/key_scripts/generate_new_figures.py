from pathlib import Path

import structlog

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.figures.fig_constants import FIGURE_DIRECTORY
from mecfs_bio.figures.figure_exporter import (
    AbstractFigureExporter,
    ValidFigureMeta,
    get_figure_destination,
)
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.key_scripts.generate_figures import FIGURE_EXPORTER

logger = structlog.get_logger()


def generate_new_figures(
    all_figure_tasks: list[Task],
    fig_dir: Path,
    exporter: AbstractFigureExporter,
):
    """
    Aim of this function is to only invoke the build system  to build figures that do not
    currently exist in the figure directory.
    """
    tasks_to_run = []
    for task in all_figure_tasks:
        meta = task.meta
        assert isinstance(meta, ValidFigureMeta)
        dst = get_figure_destination(meta=meta, fig_dir=fig_dir)
        if not dst.exists():
            tasks_to_run.append(task)

    logger.debug(
        f"Found new figure tasks: {[task.meta.asset_id for task in tasks_to_run]}"
    )
    if len(tasks_to_run) > 0:
        exporter.export(tasks_to_run, fig_dir=fig_dir)


if __name__ == "__main__":
    generate_new_figures(
        all_figure_tasks=ALL_FIGURE_TASKS,
        fig_dir=FIGURE_DIRECTORY,
        exporter=FIGURE_EXPORTER,
    )
