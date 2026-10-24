"""
Benchmark whether requested_connections speeds up a real GWAS Catalog download.

Downloads the Portas et al. sleep GWAS (GCST90824104, ~372MB) through RobustWFDownloader,
the production downloader, once per connection count, verifying the md5 hash each time.

Usage:
    pixi r python experiments/claude/download_connections_benchmark/benchmark_requested_connections.py \
        <scratch_dir> 2>&1 | tee experiments/claude/download_connections_benchmark/benchmark.log
"""

import sys
import time
from pathlib import Path

from mecfs_bio.build_system.wf.wf_downloader import RobustWFDownloader

URL = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90824001-GCST90825000/GCST90824104/harmonised/GCST90824104.h.tsv.gz"
MD5 = "8aad7906354bf7c50034d864e5d44a73"
CONNECTION_COUNTS = [8, 1]


def time_download(scratch_dir: Path, connections: int) -> float:
    dest = scratch_dir / f"connections_{connections}" / "GCST90824104.h.tsv.gz"
    assert not dest.exists(), f"{dest} already exists; use a fresh scratch dir"
    start = time.monotonic()
    RobustWFDownloader().download(
        url=URL, local_path=dest, md5_hash=MD5, requested_connections=connections
    )
    elapsed = time.monotonic() - start
    size_mb = dest.stat().st_size / 1e6
    print(
        f"RESULT connections={connections} seconds={elapsed:.1f} "
        f"size_mb={size_mb:.1f} mb_per_s={size_mb / elapsed:.2f}",
        flush=True,
    )
    dest.unlink()
    return elapsed


def main() -> None:
    scratch_dir = Path(sys.argv[1])
    for connections in CONNECTION_COUNTS:
        time_download(scratch_dir=scratch_dir, connections=connections)


if __name__ == "__main__":
    main()
