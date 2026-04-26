"""
Make the local figure directory match the manifest committed to git.

For each entry ``path -> sha256`` in the manifest, the corresponding blob is
fetched from the GitHub release (asset name = sha256) only if the local file
is missing or has a different hash. Files under the figure directory that are
not listed in the manifest are left alone unless ``prune=True`` is passed.
"""

from pathlib import Path

import structlog

from mecfs_bio.constants.gh_constants import GH_REPO_NAME
from mecfs_bio.figures.fig_constants import (
    FIGURE_DIRECTORY,
    FIGURE_GITHUB_RELEASE_TAG,
    FIGURE_MANIFEST_PATH,
)
from mecfs_bio.figures.manifest import FigureManifest, sha256_of_file
from mecfs_bio.util.github_commands.upload_download import (
    download_release_asset,
)

logger = structlog.get_logger()


def pull_figures(
    tag: str = FIGURE_GITHUB_RELEASE_TAG,
    repo_name: str = GH_REPO_NAME,
    fig_dir: Path = FIGURE_DIRECTORY,
    manifest_path: Path = FIGURE_MANIFEST_PATH,
    use_gh_cli: bool = True,
    prune: bool = False,
):
    """
    Sync the local figure directory with the manifest by downloading any
    missing or out-of-date blobs from the GitHub release.

    If ``prune`` is True, files under ``fig_dir`` that are not listed in the
    manifest are deleted.
    """
    fig_dir.mkdir(parents=True, exist_ok=True)
    manifest = FigureManifest.load(manifest_path)

    if not manifest.figures:
        logger.debug(f"Manifest {manifest_path} is empty; nothing to download.")
    for rel_path, expected_hash in manifest.figures.items():
        dest = fig_dir / rel_path
        if dest.is_file() and sha256_of_file(dest) == expected_hash:
            logger.debug(f"{rel_path} is up to date; skipping download.")
            continue
        logger.debug(f"Fetching {rel_path} (asset {expected_hash}) from release {tag}.")
        download_release_asset(
            release_tag=tag,
            repo_name=repo_name,
            asset_name=expected_hash,
            dest=dest,
            use_gh_cli=use_gh_cli,
        )
        actual = sha256_of_file(dest)
        assert actual == expected_hash, (
            f"Downloaded blob for {rel_path} hashed to {actual}, "
            f"expected {expected_hash}"
        )

    if prune:
        _prune_unmanifested(fig_dir=fig_dir, manifest=manifest)


def _prune_unmanifested(fig_dir: Path, manifest: FigureManifest) -> None:
    managed = set(manifest.figures.keys())
    for path in fig_dir.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(fig_dir).as_posix()
        if rel not in managed:
            logger.debug(f"Pruning {rel} (not in manifest).")
            path.unlink()


if __name__ == "__main__":
    pull_figures()
