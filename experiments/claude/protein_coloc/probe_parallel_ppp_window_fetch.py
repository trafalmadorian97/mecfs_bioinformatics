"""
Would fetching several proteins' dense windows at once speed up PPP colocalization?

Runs PppWindowSource.fetch_protein (production downloader, ranged reads) over batches of 12 proteins
with a thread pool of 1, 2, 4 or 8 workers. Each protein gets one 1 Mb window, on a chromosome
from a fixed mix of sizes, so every fetch costs what a typical coloc fetch costs: URL lookup,
tar header walk, then reading and decompressing one chromosome file. Every batch uses
different proteins, so no batch benefits from objects another batch warmed, and each worker
count runs twice in alternating order because repeat runs are noisy.

Reports per batch: wall time, speedup over the mean one-worker wall time, process CPU time
divided by wall time (near 1 or above means CPU-bound, so threads help less), and peak RSS.

Run: pixi r python experiments/claude/protein_coloc/probe_parallel_ppp_window_fetch.py 2>&1 | tee experiments/claude/protein_coloc/probe_parallel_ppp_window_fetch.log
"""

import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from attrs import frozen

from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval
from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
    load_ppp_protein_files,
)
from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    PppProteinFile,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_source import (
    PppWindowSource,
)
from mecfs_bio.build_system.wf.base_wf import WF, make_wf

_WORKER_ORDER = (1, 4, 2, 8, 8, 2, 4, 1)
_BATCH_SIZE = 12
_FIRST_PROTEIN = 10
# A mix of large, medium and small chromosomes (23 = X).
_CHROMS = (1, 2, 6, 11, 19, 22, 3, 8, 12, 16, 20, 23)
_WINDOW_START = 20_000_000
_WINDOW_LENGTH = 1_000_000
_RSS_POLL_SECONDS = 0.05


def _rss_mb() -> float:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) / 1024
    raise AssertionError("VmRSS missing from /proc/self/status")


class _PeakRss:
    def __init__(self) -> None:
        self.peak_mb = _rss_mb()
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._poll, daemon=True)

    def _poll(self) -> None:
        while not self._stop.wait(_RSS_POLL_SECONDS):
            self.peak_mb = max(self.peak_mb, _rss_mb())

    def __enter__(self) -> "_PeakRss":
        self._thread.start()
        return self

    def __exit__(self, *exc: object) -> None:
        self._stop.set()
        self._thread.join()


@frozen(slots=True)
class _BatchResult:
    n_workers: int
    wall_seconds: float
    cpu_seconds: float
    peak_rss_mb: float
    n_rows: int


def _fetch_one(
    source: PppWindowSource, protein: PppProteinFile, chrom: int, wf: WF
) -> int:
    window = GenomicInterval(chrom=chrom,
        start=_WINDOW_START,
        end=_WINDOW_START + _WINDOW_LENGTH,
    )
    with tempfile.TemporaryDirectory() as scratch:
        frames = source.fetch_protein(protein.oid, [window], Path(scratch), wf)
    return frames[window].height


def _run_batch(
    proteins: tuple[PppProteinFile, ...], n_workers: int, wf: WF
) -> _BatchResult:
    source = PppWindowSource(proteins=proteins)
    wall_start, cpu_start = time.monotonic(), time.process_time()
    with _PeakRss() as rss, ThreadPoolExecutor(max_workers=n_workers) as pool:
        heights = list(
            pool.map(lambda pc: _fetch_one(source, pc[0], pc[1], wf), zip(proteins, _CHROMS))
        )
    return _BatchResult(
        n_workers=n_workers,
        wall_seconds=time.monotonic() - wall_start,
        cpu_seconds=time.process_time() - cpu_start,
        peak_rss_mb=rss.peak_mb,
        n_rows=sum(heights),
    )


def main() -> None:
    proteins = load_ppp_protein_files(EUR_DISCOVERY_PPP_MANIFEST_PATH)
    wf = make_wf()
    # Log in and warm the session before timing anything.
    wf.open_synapse_range_reader(proteins[0].synid)
    print(f"baseline RSS {_rss_mb():.0f} MB")
    results: list[_BatchResult] = []
    for batch_index, n_workers in enumerate(_WORKER_ORDER):
        first = _FIRST_PROTEIN + batch_index * _BATCH_SIZE
        batch = tuple(proteins[first : first + _BATCH_SIZE])
        result = _run_batch(batch, n_workers, wf)
        results.append(result)
        print(
            f"batch {batch_index}: {n_workers} worker(s), {result.wall_seconds:.1f}s wall, "
            f"{result.wall_seconds / _BATCH_SIZE:.2f}s/protein, "
            f"cpu/wall {result.cpu_seconds / result.wall_seconds:.2f}, "
            f"peak RSS {result.peak_rss_mb:.0f} MB, {result.n_rows} rows",
            flush=True,
        )
    serial = [r.wall_seconds for r in results if r.n_workers == 1]
    serial_mean = sum(serial) / len(serial)
    print("\nsummary (mean over the two batches per worker count)")
    for n_workers in sorted(set(_WORKER_ORDER)):
        walls = [r.wall_seconds for r in results if r.n_workers == n_workers]
        mean = sum(walls) / len(walls)
        print(
            f"  {n_workers} worker(s): {mean:.1f}s per {_BATCH_SIZE} proteins, "
            f"{mean / _BATCH_SIZE:.2f}s/protein, speedup {serial_mean / mean:.2f}x"
        )


if __name__ == "__main__":
    main()
