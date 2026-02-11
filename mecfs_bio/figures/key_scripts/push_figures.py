from pathlib import Path

import structlog

"""
Merge the contents of the local figure directory with figures stored on Github, then upload to Github the result of this merge.
"""
from mecfs_bio.constants.gh_constants import GH_REPO_NAME

from mecfs_bio.figures.fig_constants import (
    FIGURE_DIRECTORY,
    FIGURE_GITHUB_RELEASE_TAG,
    FIGURES_ARCHIVE_TITLE,
)
from mecfs_bio.figures.key_scripts.pull_figures import pull_figures
from mecfs_bio.util.github_commands.upload_download import (
    release_from_dir,
)

logger = structlog.get_logger()


def push_figures(
    tag: str = FIGURE_GITHUB_RELEASE_TAG,
    repo_name: str = GH_REPO_NAME,
    fig_dir: Path = FIGURE_DIRECTORY,
):
    """
    Merge the contents of the local figure directory with figures stored on Github, then upload to Github the result of this merge.
    """
    fig_dir.mkdir(parents=True, exist_ok=True)
    pull_figures(
        tag=tag,
        repo_name=repo_name,
        fig_dir=fig_dir,
    )
    release_from_dir(
        release_tag=tag,
        dir_path=fig_dir,
        title=FIGURES_ARCHIVE_TITLE,
        repo_name=repo_name,
        draft=False,
    )
    logger.debug("Figures successfully uploaded")


if __name__ == "__main__":
    push_figures()
