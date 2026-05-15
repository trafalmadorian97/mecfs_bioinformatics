"""
Update the figure manifest to reflect the local figure directory and upload any
new blobs to the GitHub release.

The local figure directory is treated as the source of truth for adds and
updates: every file under it is hashed and recorded in the manifest. Blobs are
content-addressed --- each unique hash is uploaded at most once and never
overwritten, so updates by different collaborators do not clobber each other
on the release. Manifest changes are written to disk for the user to commit.

By default, entries already present in the manifest but absent from the local
figure directory are left in place (additive behaviour). Pass ``prune=True``
to drop those entries from the manifest --- the upstream blobs remain on the
release for any past commit that still references them.
"""

from pathlib import Path
from typing import Sequence

import structlog

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.constants.gh_constants import GH_REPO_NAME
from mecfs_bio.figures.fig_constants import (
    FIGURE_DIRECTORY,
    FIGURE_GITHUB_RELEASE_TAG,
    FIGURE_MANIFEST_PATH,
    FIGURES_ARCHIVE_TITLE,
)
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.manifest import (
    FigureManifest,
    scan_figure_dir,
)
from mecfs_bio.figures.manifest_validation import validate_manifest_subset_of_tasks
from mecfs_bio.util.github_commands.upload_download import (
    list_release_asset_names,
    upload_blob_to_release,
)

logger = structlog.get_logger()


def push_figures(
    tag: str = FIGURE_GITHUB_RELEASE_TAG,
    repo_name: str = GH_REPO_NAME,
    fig_dir: Path = FIGURE_DIRECTORY,
    manifest_path: Path = FIGURE_MANIFEST_PATH,
    title: str = FIGURES_ARCHIVE_TITLE,
    prune: bool = False,
    figure_tasks: Sequence[Task] = ALL_FIGURE_TASKS,
):
    """
    Update the manifest from the local figure directory and upload any new
    blobs to the GitHub release.
    """
    fig_dir.mkdir(parents=True, exist_ok=True)

    old_manifest = FigureManifest.load(manifest_path)
    local_manifest = scan_figure_dir(fig_dir)

    new_manifest = _merge_manifests(old=old_manifest, local=local_manifest, prune=prune)

    # Fail fast if the manifest we are about to write references files that
    # no task in figure_tasks produces --- those would silently rot.
    validate_manifest_subset_of_tasks(
        manifest=new_manifest, tasks=figure_tasks, fig_dir=fig_dir
    )

    remote_assets = list_release_asset_names(release_tag=tag, repo_name=repo_name)

    new_hashes = new_manifest.hashes() - remote_assets
    if not new_hashes:
        logger.debug("No new blobs to upload to release.")
    for rel_path, sha in local_manifest.figures.items():
        if sha not in new_hashes:
            continue
        src = fig_dir / rel_path
        logger.debug(f"Uploading blob {sha} (from {rel_path}) to release {tag}.")
        upload_blob_to_release(
            release_tag=tag,
            repo_name=repo_name,
            asset_name=sha,
            src_path=src,
            title=title,
        )
        new_hashes.discard(sha)

    if new_hashes:
        # A hash listed in the new manifest is missing both from the release
        # and from the local figure directory --- this can only happen if
        # the manifest already referenced a blob the user does not have a
        # local copy of and which was never uploaded.
        raise RuntimeError(
            f"Manifest references hashes that are neither on the release nor "
            f"available locally: {sorted(new_hashes)}"
        )

    new_manifest.save(manifest_path)
    logger.debug(f"Manifest written to {manifest_path}.")
    logger.debug("Push complete. Commit the manifest to record the change.")


def _merge_manifests(
    old: FigureManifest, local: FigureManifest, prune: bool
) -> FigureManifest:
    """
    Build the manifest to write back: local entries always win for paths
    present locally; paths only in ``old`` are dropped if ``prune`` is set,
    otherwise they are preserved.
    """
    merged: dict[str, str] = {} if prune else dict(old.figures)
    merged.update(local.figures)
    return FigureManifest(figures=merged)


if __name__ == "__main__":
    push_figures()
