"""
Make the local figure directory match the manifest committed to git.

For each entry ``path -> sha256`` in the manifest, the corresponding blob is
fetched from the GitHub release (asset name = sha256) only if the local file
is missing or has a different hash. Files under the figure directory that are
not listed in the manifest are left alone unless ``prune=True`` is passed.
"""

import shutil
import tempfile
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
    download_release_assets_chunked,
)

logger = structlog.get_logger()

# Chunked-download tuning. Each chunk is a single `gh release download` call
# with one repeated --pattern per asset; a small pool of chunks runs
# concurrently. This amortizes gh's per-call startup + asset-listing overhead
# (benchmarked ~6x faster than per-asset parallel downloads, with far fewer
# API calls); see experiments/claude/figure_download_benchmark.
DEFAULT_CHUNK_SIZE = 20
DEFAULT_CHUNK_WORKERS = 4


def pull_figures(
    tag: str = FIGURE_GITHUB_RELEASE_TAG,
    repo_name: str = GH_REPO_NAME,
    fig_dir: Path = FIGURE_DIRECTORY,
    manifest_path: Path = FIGURE_MANIFEST_PATH,
    use_gh_cli: bool = True,
    prune: bool = False,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    chunk_workers: int = DEFAULT_CHUNK_WORKERS,
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

    to_download: list[tuple[Path, str]] = []
    for rel_path, expected_hash in manifest.figures.items():
        dest = fig_dir / rel_path
        if dest.is_file() and sha256_of_file(dest) == expected_hash:
            logger.debug(f"{rel_path} is up to date; skipping download.")
            continue
        to_download.append((rel_path, expected_hash))

    if to_download:
        _download_blobs(
            to_download=to_download,
            tag=tag,
            repo_name=repo_name,
            fig_dir=fig_dir,
            use_gh_cli=use_gh_cli,
            chunk_size=chunk_size,
            chunk_workers=chunk_workers,
        )

    if prune:
        _prune_unmanifested(fig_dir=fig_dir, manifest=manifest)


def _download_blobs(
    to_download: list[tuple[Path, str]],
    tag: str,
    repo_name: str,
    fig_dir: Path,
    use_gh_cli: bool,
    chunk_size: int,
    chunk_workers: int,
) -> None:
    """
    Stage every needed blob (deduplicated by hash) into a temp directory, then
    place each into ``fig_dir`` at its manifest path, verifying its hash.

    Blobs are keyed by content hash, so two manifest entries that share a hash
    are downloaded once and copied to both destinations.
    """
    # dict.fromkeys deduplicates while preserving order.
    unique_hashes = list(dict.fromkeys(h for _, h in to_download))
    logger.debug(
        f"Downloading {len(unique_hashes)} unique blobs for "
        f"{len(to_download)} figure(s)."
    )
    with tempfile.TemporaryDirectory() as staging:
        staging_dir = Path(staging)
        if use_gh_cli:
            download_release_assets_chunked(
                release_tag=tag,
                repo_name=repo_name,
                asset_names=unique_hashes,
                dest_dir=staging_dir,
                chunk_size=chunk_size,
                chunk_workers=chunk_workers,
            )
        else:
            _stage_assets_individually(
                asset_names=unique_hashes,
                staging_dir=staging_dir,
                tag=tag,
                repo_name=repo_name,
            )
        _place_staged_blobs(
            to_download=to_download, staging_dir=staging_dir, fig_dir=fig_dir
        )


def _stage_assets_individually(
    asset_names: list[str],
    staging_dir: Path,
    tag: str,
    repo_name: str,
) -> None:
    """
    Fallback path for ``use_gh_cli=False``: download each blob by URL (no gh,
    no auth) into the staging dir. The chunked --pattern trick is gh-only.
    """
    for asset_name in asset_names:
        download_release_asset(
            release_tag=tag,
            repo_name=repo_name,
            asset_name=asset_name,
            dest=staging_dir / asset_name,
            use_gh_cli=False,
        )


def _place_staged_blobs(
    to_download: list[tuple[Path, str]],
    staging_dir: Path,
    fig_dir: Path,
) -> None:
    for rel_path, expected_hash in to_download:
        src = staging_dir / expected_hash
        dest = fig_dir / rel_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        # copyfile copies data only (no chmod/utime), avoiding cross-device and
        # DrvFs metadata pitfalls.
        shutil.copyfile(src, dest)
        actual = sha256_of_file(dest)
        assert actual == expected_hash, (
            f"Blob for {rel_path} hashed to {actual}, expected {expected_hash}"
        )


def _prune_unmanifested(fig_dir: Path, manifest: FigureManifest) -> None:
    managed = set(manifest.figures.keys())
    for path in fig_dir.rglob("*"):
        if not path.is_file():
            continue
        rel = path.relative_to(fig_dir)
        if rel not in managed:
            logger.debug(f"Pruning {rel} (not in manifest).")
            path.unlink()


if __name__ == "__main__":
    pull_figures()
