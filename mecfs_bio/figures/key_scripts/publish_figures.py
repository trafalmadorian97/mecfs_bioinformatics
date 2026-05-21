"""
End-to-end "the figure system is out of date" workflow.

Runs the three steps a contributor wants in sequence:
  1. ``generate_new_figures`` --- invoke the build system for any figure
     task in ALL_FIGURE_TASKS whose output is not yet on disk, and copy the
     results into the figure directory.
  2. ``prune_orphan_figures`` --- drop manifest entries (and local files)
     for figures no task in ALL_FIGURE_TASKS produces. If any such entry
     is still referenced in the documentation, the script raises before
     touching anything; the user resolves the conflict by restoring the
     task or removing the doc reference.
  3. ``push_figures`` --- rehash the figure directory, update the
     committed manifest, and upload any new content-addressed blobs to
     the GitHub release.

After this script finishes, commit the updated
``mecfs_bio/figures/figures_manifest.json`` to record the change.
"""

import structlog

from mecfs_bio.figures.fig_constants import (
    DOCS_DIRECTORY,
    FIGURE_DIRECTORY,
    FIGURE_MANIFEST_PATH,
)
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.key_scripts.generate_figures import FIGURE_EXPORTER
from mecfs_bio.figures.key_scripts.generate_new_figures import generate_new_figures
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.orphan_pruning import prune_orphan_figures

logger = structlog.get_logger()


def publish_figures() -> None:
    logger.info("Generating any missing figures from ALL_FIGURE_TASKS.")
    generate_new_figures(
        all_figure_tasks=ALL_FIGURE_TASKS,
        fig_dir=FIGURE_DIRECTORY,
        exporter=FIGURE_EXPORTER,
    )
    logger.info("Pruning orphan figures (entries not produced by ALL_FIGURE_TASKS).")
    prune_orphan_figures(
        manifest_path=FIGURE_MANIFEST_PATH,
        docs_dir=DOCS_DIRECTORY,
        fig_dir=FIGURE_DIRECTORY,
        figure_tasks=ALL_FIGURE_TASKS,
    )
    logger.info("Updating manifest and pushing new blobs to the GitHub release.")
    push_figures(figure_tasks=ALL_FIGURE_TASKS)
    logger.info(
        "Done. Commit mecfs_bio/figures/figures_manifest.json to record the change."
    )


if __name__ == "__main__":
    publish_figures()
