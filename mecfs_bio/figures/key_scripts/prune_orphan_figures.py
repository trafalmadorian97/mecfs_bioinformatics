"""
Remove figure-manifest entries that no task in ALL_FIGURE_TASKS produces.

Refuses to prune when the docs still reference the orphan. Run as a
standalone tool, or rely on ``publish_new_figures`` which calls the same
function between ``generate_new_figures`` and ``push_figures``.
"""

from mecfs_bio.figures.fig_constants import (
    DOCS_DIRECTORY,
    FIGURE_DIRECTORY,
    FIGURE_MANIFEST_PATH,
)
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.orphan_pruning import prune_orphan_figures


def main() -> None:
    prune_orphan_figures(
        manifest_path=FIGURE_MANIFEST_PATH,
        docs_dir=DOCS_DIRECTORY,
        fig_dir=FIGURE_DIRECTORY,
        figure_tasks=ALL_FIGURE_TASKS,
    )


if __name__ == "__main__":
    main()
