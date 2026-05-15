import json
import shutil
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
            for file_path in dir_path.rglob("*"):
                archive_path = file_path.relative_to(dir_path)
                zipf.write(file_path, arcname=archive_path)
                logger.debug(f"Added {archive_path} to figures release")

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

        tmp_output_path = Path(tmpdirname) / "tmp_dir"
        tmp_output_path.mkdir(parents=True, exist_ok=True)
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
            zip_ref.extractall(tmp_output_path)
            shutil.copytree(tmp_output_path, dir_path, dirs_exist_ok=True)


def download_release_to_dir_no_auth(
    release_tag: str, dir_path: Path, repo_name: str, title: str
) -> None:
    """
    As above, but don't use the Github command line tool, so do not require authorization.
    """
    assert dir_path.is_dir()
    with tempfile.TemporaryDirectory() as tmpdirname:
        zip_path = Path(tmpdirname) / f"{release_tag}.zip"
        tmp_output_path = Path(tmpdirname) / "tmp_dir"
        tmp_output_path.mkdir(parents=True, exist_ok=True)
        robust_download_with_aria(
            md5sum=None,
            dest=zip_path,
            url=get_release_url(
                release_tag=release_tag, repo_name=repo_name, title=title
            ),
        )

        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(tmp_output_path)
            shutil.copytree(tmp_output_path, dir_path, dirs_exist_ok=True)


def get_release_url(release_tag: str, repo_name: str, title: str) -> str:
    return f"https://github.com/{repo_name}/releases/download/{release_tag}/{title}.zip"


def delete_release(release_tag: str, repo_name: str, cleanup_tag: bool = True) -> None:
    """
    Delete a release. No-op if the release does not exist. By default, the
    underlying git tag is also removed so the release name can be reused.
    """
    if not does_release_exist(repo_name=repo_name, release_tag=release_tag):
        return
    cmd = [
        "gh",
        "release",
        "delete",
        release_tag,
        "-R",
        repo_name,
        "--yes",
    ]
    if cleanup_tag:
        cmd.append("--cleanup-tag")
    execute_command(cmd)


def list_release_asset_names(release_tag: str, repo_name: str) -> set[str]:
    """
    Return the set of asset names attached to a release. Returns an empty set
    if the release does not exist.
    """
    if not does_release_exist(repo_name=repo_name, release_tag=release_tag):
        return set()
    cmd = [
        "gh",
        "release",
        "view",
        release_tag,
        "-R",
        repo_name,
        "--json",
        "assets",
    ]
    output = execute_command(cmd)
    payload = json.loads(output)
    return {asset["name"] for asset in payload.get("assets", [])}


def ensure_release_exists(
    release_tag: str, repo_name: str, title: str | None = None
) -> None:
    """
    Idempotently create the GitHub release if it does not already exist.

    Callers that intend to upload many blobs should invoke this once up
    front; subsequent ``upload_blob_to_release`` calls can then skip the
    per-blob existence check (which previously doubled the number of
    ``gh`` invocations).
    """
    if does_release_exist(release_tag=release_tag, repo_name=repo_name):
        return
    cmd = [
        "gh",
        "release",
        "create",
        release_tag,
        "--title",
        title if title is not None else release_tag,
        "-R",
        repo_name,
    ]
    execute_command(cmd=cmd)


def upload_blob_to_release(
    release_tag: str,
    repo_name: str,
    asset_name: str,
    src_path: Path,
) -> None:
    """
    Upload ``src_path`` to the release as an asset named ``asset_name``.

    GitHub uses the local filename as the asset name, so the file is staged
    under ``asset_name`` in a temp dir before upload. The release must
    already exist; call ``ensure_release_exists`` once before the first
    upload.
    """
    assert src_path.is_file()
    with tempfile.TemporaryDirectory() as tmpdir:
        staged = Path(tmpdir) / asset_name
        shutil.copy(src_path, staged)
        cmd = [
            "gh",
            "release",
            "upload",
            release_tag,
            str(staged),
            "-R",
            repo_name,
        ]
        execute_command(cmd=cmd)


def download_release_asset(
    release_tag: str,
    repo_name: str,
    asset_name: str,
    dest: Path,
    use_gh_cli: bool = True,
) -> None:
    """
    Download a single named asset from a release to ``dest``.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if use_gh_cli:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_dir = Path(tmpdir)
            cmd = [
                "gh",
                "release",
                "download",
                release_tag,
                "--pattern",
                asset_name,
                "--dir",
                str(tmp_dir),
                "-R",
                repo_name,
            ]
            execute_command(cmd)
            downloaded = tmp_dir / asset_name
            assert downloaded.is_file(), (
                f"gh did not produce {downloaded} for asset {asset_name}"
            )
            shutil.move(str(downloaded), dest)
    else:
        url = (
            f"https://github.com/{repo_name}/releases/download/"
            f"{release_tag}/{asset_name}"
        )
        robust_download_with_aria(md5sum=None, dest=dest, url=url)
