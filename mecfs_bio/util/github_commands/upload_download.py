import tempfile
import zipfile
from pathlib import Path
from subprocess import CalledProcessError

import structlog

from mecfs_bio.util.download.robust_download import robust_download_with_aria
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()


def release_from_dir(
    release_tag: str, dir_path: Path, title: str, draft: bool, repo_name: str
):
    """
    Create a github release by zipping the contents of a directory
    """
    assert dir_path.is_dir()

    with tempfile.TemporaryDirectory() as tmpdirname:
        temp_path = Path(tmpdirname)
        zip_path = temp_path / f"{title}.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
            for file_path in dir_path.glob("*"):
                zipf.write(file_path, arcname=file_path.name)

        if not does_release_exist(release_tag=release_tag, repo_name=repo_name):
            cmd = [
                "gh",
                "release",
                "create",
                release_tag,
                str(zip_path),
                "--title",
                title,
                "-R",
                repo_name,
            ]
            if draft:
                cmd.append("--draft")
        else:
            cmd = [
                "gh",
                "release",
                "upload",
                "--clobber",
                release_tag,
                str(zip_path),
            ]
        execute_command(cmd=cmd)


def does_release_exist(repo_name: str, release_tag: str):
    cmd = ["gh", "release", "view", "-R", repo_name, release_tag]
    try:
        execute_command(cmd)
        return True
    except CalledProcessError as e:
        logger.debug(f"Got status code {e.returncode}")
        if (
            e.returncode == 1
        ):  # Return codes are ambiguous, but 1 typically means release not found: https://cli.github.com/manual/gh_help_exit-codes#:~:text=If%20a%20command%20completes%20successfully,exit%20code%20will%20be%204
            return False
        raise e


def download_release_to_dir(release_tag: str, dir_path: Path, repo_name: str):
    """
    Download a github release, which is assumed to be a zip file. Extract it into a directory.
    """
    assert dir_path.is_dir()
    with tempfile.TemporaryDirectory() as tmpdirname:
        zip_path = Path(tmpdirname) / f"{release_tag}.zip"
        cmd = [
            "gh",
            "release",
            "download",
            release_tag,
            "--output",
            str(zip_path),
            "-R",
            repo_name,
        ]
        execute_command(cmd)
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(dir_path)


def download_release_to_dir_no_auth(
    release_tag: str, dir_path: Path, repo_name: str, title: str
) -> None:
    """
    As above, but don't use the Github command line tool, so do not require authorization.
    """
    assert dir_path.is_dir()
    with tempfile.TemporaryDirectory() as tmpdirname:
        zip_path = Path(tmpdirname) / f"{release_tag}.zip"
        robust_download_with_aria(
            md5sum=None,
            dest=zip_path,
            url=get_release_url(
                release_tag=release_tag, repo_name=repo_name, title=title
            ),
        )

        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(dir_path)


def get_release_url(release_tag: str, repo_name: str, title: str) -> str:
    return f"https://github.com/{repo_name}/releases/download/{release_tag}/{title}.zip"
