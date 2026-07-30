"""
Benchmark two strategies for pulling the figure blobs from the GitHub release.

Strategy A (current production path): one `gh release download --pattern <hash>`
invocation per asset, run in a thread pool. Each invocation re-lists the
release's assets before downloading its one blob.

Strategy B (chunked batch): split the manifest's assets into chunks and issue a
single `gh release download` per chunk with one repeated `-p <hash>` per asset,
so a chunk is located with a single asset-list lookup. Chunks themselves are run
in a small thread pool.

Both strategies download exactly the manifest's assets (never the release's
orphaned/old blobs), so the manifest-consistency guarantee is preserved.

Each configuration is timed into a fresh temp dir and verified to have produced
every requested blob. Wall-clock is reported as the median over repetitions to
blunt network variance.

Run: pixi r python experiments/claude/figure_download_benchmark/benchmark.py
"""

import statistics
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import structlog
from attrs import frozen

from mecfs_bio.constants.gh_constants import GH_REPO_NAME
from mecfs_bio.figures.fig_constants import (
    FIGURE_GITHUB_RELEASE_TAG,
    FIGURE_MANIFEST_PATH,
)
from mecfs_bio.figures.manifest import FigureManifest
from mecfs_bio.util.github_commands.upload_download import download_release_asset
from mecfs_bio.util.subproc.run_command import execute_command_with_retries

logger = structlog.get_logger()

REPETITIONS = 3


def _unique_assets() -> list[str]:
    manifest = FigureManifest.load(FIGURE_MANIFEST_PATH)
    # Preserve order but de-duplicate (two figures can share a hash).
    seen: dict[str, None] = {}
    for h in manifest.figures.values():
        seen.setdefault(h, None)
    return list(seen.keys())


def download_per_asset_parallel(
    assets: list[str], dest: Path, workers: int
) -> None:
    def _one(asset: str) -> None:
        download_release_asset(
            release_tag=FIGURE_GITHUB_RELEASE_TAG,
            repo_name=GH_REPO_NAME,
            asset_name=asset,
            dest=dest / asset,
            use_gh_cli=True,
        )

    with ThreadPoolExecutor(max_workers=min(workers, len(assets))) as pool:
        list(pool.map(_one, assets))


def download_chunked_batch(
    assets: list[str], dest: Path, chunk_size: int, chunk_workers: int
) -> None:
    chunks = [assets[i : i + chunk_size] for i in range(0, len(assets), chunk_size)]

    def _one_chunk(chunk: list[str]) -> None:
        cmd = ["gh", "release", "download", FIGURE_GITHUB_RELEASE_TAG]
        for asset in chunk:
            cmd += ["--pattern", asset]
        cmd += ["--dir", str(dest), "--clobber", "-R", GH_REPO_NAME]
        execute_command_with_retries(cmd)

    with ThreadPoolExecutor(max_workers=min(chunk_workers, len(chunks))) as pool:
        list(pool.map(_one_chunk, chunks))


@frozen
class Config:
    name: str
    run: object  # callable (assets, dest) -> None

    def execute(self, assets: list[str], dest: Path) -> None:
        self.run(assets, dest)  # type: ignore[operator]


def _verify(dest: Path, assets: list[str]) -> None:
    missing = [a for a in assets if not (dest / a).is_file()]
    assert not missing, f"{len(missing)} assets missing after download: {missing[:3]}"


def time_config(config: Config, assets: list[str]) -> float:
    with tempfile.TemporaryDirectory() as tmp:
        dest = Path(tmp)
        start = time.monotonic()
        config.execute(assets, dest)
        elapsed = time.monotonic() - start
        _verify(dest, assets)
    return elapsed


def main() -> None:
    assets = _unique_assets()
    logger.info(f"Benchmarking with {len(assets)} unique assets, {REPETITIONS} reps.")

    configs = [
        Config(
            "per_asset_parallel_w8",
            lambda a, d: download_per_asset_parallel(a, d, workers=8),
        ),
        Config(
            "chunked_c20_w4",
            lambda a, d: download_chunked_batch(a, d, chunk_size=20, chunk_workers=4),
        ),
        Config(
            "chunked_c20_w1",
            lambda a, d: download_chunked_batch(a, d, chunk_size=20, chunk_workers=1),
        ),
        Config(
            "chunked_all_w1_bulk",
            lambda a, d: download_chunked_batch(
                a, d, chunk_size=len(a), chunk_workers=1
            ),
        ),
    ]

    results: dict[str, list[float]] = {}
    for config in configs:
        times = []
        for rep in range(REPETITIONS):
            elapsed = time_config(config, assets)
            logger.info(f"{config.name} rep {rep + 1}: {elapsed:.1f}s")
            times.append(elapsed)
        results[config.name] = times

    print("\n=== SUMMARY (wall-clock seconds) ===")
    print(f"{'config':<26} {'median':>8} {'min':>8} {'max':>8}  runs")
    for name, times in results.items():
        print(
            f"{name:<26} {statistics.median(times):>8.1f} "
            f"{min(times):>8.1f} {max(times):>8.1f}  "
            f"{[round(t, 1) for t in times]}"
        )


if __name__ == "__main__":
    main()
