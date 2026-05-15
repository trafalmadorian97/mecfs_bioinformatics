"""
End-to-end "I added a new figure" workflow.

Runs the two steps a contributor wants in sequence:
  1. ``generate_new_figures`` --- invoke the build system for any figure
     task in ALL_FIGURE_TASKS whose output is not yet on disk, and copy the
     results into the figure directory.
  2. ``push_figures`` --- rehash the figure directory, update the committed
     manifest, and upload any new content-addressed blobs to the GitHub
     release.

After this script finishes, commit the updated
``mecfs_bio/figures/figures_manifest.json`` to record the change.
"""

import structlog

from mecfs_bio.figures.fig_constants import FIGURE_DIRECTORY
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.key_scripts.generate_figures import FIGURE_EXPORTER
from mecfs_bio.figures.key_scripts.generate_new_figures import generate_new_figures
from mecfs_bio.figures.key_scripts.push_figures import push_figures

logger = structlog.get_logger()


def publish_new_figures() -> None:
    logger.info("Generating any missing figures from ALL_FIGURE_TASKS.")
    generate_new_figures(
        all_figure_tasks=ALL_FIGURE_TASKS,
        fig_dir=FIGURE_DIRECTORY,
        exporter=FIGURE_EXPORTER,
    )
    logger.info("Updating manifest and pushing new blobs to the GitHub release.")
    push_figures(figure_tasks=ALL_FIGURE_TASKS)
    logger.info(
        "Done. Commit mecfs_bio/figures/figures_manifest.json to record the change."
    )


if __name__ == "__main__":
    publish_new_figures()
