"""
End-to-end smoke test for the manifest-based figure system.

Exercises real ``gh release`` upload/download against a throwaway tag.
Requires ``gh`` to be installed and authenticated, and write access to
``GH_REPO_NAME``. The release (and its underlying git tag) is deleted at the
start and end of the run, so the script can be invoked repeatedly.

Run via ``pixi r python experiments/tralfamadorian97/smoke_test_figure_manifest.py``.
"""

import shutil
import sys
import tempfile
from pathlib import Path

import structlog

from mecfs_bio.constants.gh_constants import GH_REPO_NAME
from mecfs_bio.figures.key_scripts.pull_figures import pull_figures
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.manifest import FigureManifest, sha256_of_file
from mecfs_bio.util.github_commands.upload_download import (
    delete_release,
    list_release_asset_names,
)

logger = structlog.get_logger()

SMOKE_TAG = "figures-smoke-test"
SMOKE_TITLE = "figures-smoke-test"


def _write(path: Path, contents: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(contents)


def _check(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(f"smoke test failed: {message}")
    logger.debug(f"OK: {message}")


def _push(fig_dir: Path, manifest_path: Path, prune: bool = False) -> None:
    push_figures(
        tag=SMOKE_TAG,
        repo_name=GH_REPO_NAME,
        fig_dir=fig_dir,
        manifest_path=manifest_path,
        title=SMOKE_TITLE,
        prune=prune,
    )


def _pull(fig_dir: Path, manifest_path: Path, prune: bool = False) -> None:
    pull_figures(
        tag=SMOKE_TAG,
        repo_name=GH_REPO_NAME,
        fig_dir=fig_dir,
        manifest_path=manifest_path,
        prune=prune,
    )


def run_smoke_test() -> None:
    delete_release(release_tag=SMOKE_TAG, repo_name=GH_REPO_NAME)
    try:
        with tempfile.TemporaryDirectory() as tmp:
            workspace = Path(tmp)
            fig_dir = workspace / "figs_a"
            manifest_path = workspace / "manifest_a.json"

            logger.info("== step 1: initial push of two figures ==")
            _write(fig_dir / "fig1.png", b"first figure")
            _write(fig_dir / "sub/fig2.html", b"<html>second</html>")
            _push(fig_dir=fig_dir, manifest_path=manifest_path)

            manifest = FigureManifest.load(manifest_path)
            _check(
                set(manifest.figures.keys()) == {"fig1.png", "sub/fig2.html"},
                "manifest contains both figures after initial push",
            )
            remote = list_release_asset_names(
                release_tag=SMOKE_TAG, repo_name=GH_REPO_NAME
            )
            _check(
                remote == manifest.hashes(),
                "release assets exactly match manifest hashes",
            )

            logger.info("== step 2: pull into fresh dir reproduces files ==")
            fig_dir_b = workspace / "figs_b"
            _pull(fig_dir=fig_dir_b, manifest_path=manifest_path)
            _check(
                (fig_dir_b / "fig1.png").read_bytes() == b"first figure",
                "fig1.png contents match after pull",
            )
            _check(
                (fig_dir_b / "sub/fig2.html").read_bytes()
                == b"<html>second</html>",
                "sub/fig2.html contents match after pull",
            )

            logger.info(
                "== step 3: re-push with no changes uploads no new blobs =="
            )
            assets_before = list_release_asset_names(
                release_tag=SMOKE_TAG, repo_name=GH_REPO_NAME
            )
            _push(fig_dir=fig_dir, manifest_path=manifest_path)
            assets_after = list_release_asset_names(
                release_tag=SMOKE_TAG, repo_name=GH_REPO_NAME
            )
            _check(
                assets_before == assets_after,
                "no-op push leaves release assets unchanged",
            )

            logger.info(
                "== step 4: update fig1 -- new blob added, old blob retained =="
            )
            old_fig1_hash = sha256_of_file(fig_dir / "fig1.png")
            _write(fig_dir / "fig1.png", b"first figure -- v2")
            _push(fig_dir=fig_dir, manifest_path=manifest_path)
            manifest = FigureManifest.load(manifest_path)
            new_fig1_hash = manifest.figures["fig1.png"]
            _check(
                new_fig1_hash != old_fig1_hash,
                "manifest now points to new hash for fig1.png",
            )
            remote = list_release_asset_names(
                release_tag=SMOKE_TAG, repo_name=GH_REPO_NAME
            )
            _check(
                old_fig1_hash in remote,
                "old fig1 blob retained on release (append-only)",
            )
            _check(
                new_fig1_hash in remote,
                "new fig1 blob uploaded to release",
            )

            logger.info(
                "== step 5: pull onto stale dir picks up the update =="
            )
            _check(
                (fig_dir_b / "fig1.png").read_bytes() == b"first figure",
                "stale dir still has old fig1 contents",
            )
            _pull(fig_dir=fig_dir_b, manifest_path=manifest_path)
            _check(
                (fig_dir_b / "fig1.png").read_bytes()
                == b"first figure -- v2",
                "stale dir has updated fig1 contents after pull",
            )

            logger.info(
                "== step 6: delete fig2 locally, push without prune retains it =="
            )
            (fig_dir / "sub/fig2.html").unlink()
            _push(fig_dir=fig_dir, manifest_path=manifest_path, prune=False)
            manifest = FigureManifest.load(manifest_path)
            _check(
                "sub/fig2.html" in manifest.figures,
                "non-prune push keeps fig2 in manifest",
            )

            logger.info(
                "== step 7: push with prune drops fig2 from manifest =="
            )
            _push(fig_dir=fig_dir, manifest_path=manifest_path, prune=True)
            manifest = FigureManifest.load(manifest_path)
            _check(
                "sub/fig2.html" not in manifest.figures,
                "prune push removes fig2 entry from manifest",
            )
            _check(
                "fig1.png" in manifest.figures,
                "prune push keeps fig1 entry in manifest",
            )

            logger.info(
                "== step 8: pull with prune deletes local fig2 =="
            )
            _check(
                (fig_dir_b / "sub/fig2.html").exists(),
                "fig2 still on disk in stale dir before prune-pull",
            )
            _pull(
                fig_dir=fig_dir_b, manifest_path=manifest_path, prune=True
            )
            _check(
                not (fig_dir_b / "sub/fig2.html").exists(),
                "prune-pull deletes local fig2",
            )
            _check(
                (fig_dir_b / "fig1.png").exists(),
                "prune-pull preserves local fig1",
            )

            shutil.rmtree(fig_dir_b, ignore_errors=True)

        logger.info("smoke test PASSED")
    finally:
        delete_release(release_tag=SMOKE_TAG, repo_name=GH_REPO_NAME)


if __name__ == "__main__":
    try:
        run_smoke_test()
    except AssertionError as exc:
        logger.error(str(exc))
        sys.exit(1)
