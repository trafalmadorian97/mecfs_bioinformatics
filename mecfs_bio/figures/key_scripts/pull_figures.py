import shutil
import tempfile
from pathlib import Path

import structlog

from mecfs_bio.constants.gh_constants import GH_REPO_NAME
from mecfs_bio.figures.fig_constants import FIGURE_DIRECTORY, FIGURE_GITHUB_RELEASE_TAG
from mecfs_bio.util.github_commands.upload_download import (
    does_release_exist,
    download_release_to_dir,
)

logger = structlog.get_logger()


def pull_figures(
    tag: str = FIGURE_GITHUB_RELEASE_TAG,
    repo_name: str = GH_REPO_NAME,
    fig_dir: Path = FIGURE_DIRECTORY,
):
    """
    Download figures from Github, merging them with contents of the local figure directory.
    """
    fig_dir.mkdir(parents=True, exist_ok=True)
    if not does_release_exist(repo_name=repo_name, release_tag=tag):
        logger.debug(
            f"No release with the tag '{tag}' exists in repository {repo_name}.  Nothing to download."
        )
        return
    with tempfile.TemporaryDirectory() as tmpdir:
        staging_dir = Path(tmpdir)
        logger.debug(f"Downloading {tag} to {staging_dir}")
        download_release_to_dir(
            release_tag=tag,
            dir_path=staging_dir,
            repo_name=repo_name,
        )
        logger.debug("download complete")
        logger.debug(
            "Overlaying existing figures on downloaded figures in staging directory."
        )
        shutil.copytree(fig_dir, staging_dir, dirs_exist_ok=True)
        logger.debug(
            f"Copying figures back from staging directory to figures directory {fig_dir}."
        )
        shutil.copytree(staging_dir, fig_dir, dirs_exist_ok=True)


if __name__ == "__main__":
    pull_figures()
