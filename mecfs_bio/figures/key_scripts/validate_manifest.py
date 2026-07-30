"""
Check that the figure manifest committed to git stays a subset of the
figures produced by ALL_FIGURE_TASKS.

Run as a standalone script (e.g. in CI) to catch drift between the two
sources of truth before it reaches the documentation build.
"""

import structlog

from mecfs_bio.figures.fig_constants import FIGURE_DIRECTORY, FIGURE_MANIFEST_PATH
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.manifest import FigureManifest
from mecfs_bio.figures.manifest_validation import validate_manifest_subset_of_tasks

logger = structlog.get_logger()


def main() -> None:
    manifest = FigureManifest.load(FIGURE_MANIFEST_PATH)
    validate_manifest_subset_of_tasks(
        manifest=manifest, tasks=ALL_FIGURE_TASKS, fig_dir=FIGURE_DIRECTORY
    )
    logger.info(
        f"Manifest {FIGURE_MANIFEST_PATH} is a subset of figures produced by "
        f"ALL_FIGURE_TASKS ({len(manifest.figures)} entries checked)."
    )


if __name__ == "__main__":
    main()
