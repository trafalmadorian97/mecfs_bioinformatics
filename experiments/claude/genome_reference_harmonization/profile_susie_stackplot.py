"""
Memory-profile the SUSIE stackplot's LD-heatmap path to explain the OOM seen when
rerunning the DecodeME fine-mapping (chr20 47.0-48.2M, palindromes kept).

The stackplot loads an n x n LD matrix and, because it is configured with
heatmap_bin_options=None, renders it UNBINNED: get_array_and_edges_for_ld_heatmap runs
xarray da.groupby("x").mean().groupby("y").mean(), then pcolormesh draws ~n^2 quads.
The raw matrix is small (~118 MB, n~3852), so any blow-up is in those two stages. This
script measures peak anonymous RSS around each stage in isolation, so we can see which
one is responsible without running the whole plot.

RssAnon (from /proc/self/status) excludes file-backed pages and is the same metric the
kernel OOM killer reported. Each stage runs and prints the peak RSS reached during it;
if a stage itself OOMs, the log stops at that stage's banner, which is the answer.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.profile_susie_stackplot \
    2>&1 | tee experiments/claude/genome_reference_harmonization/profile_susie_stackplot.log
"""

import gc
import threading
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    get_array_and_edges_for_ld_heatmap,
)

SUSIE_DIR = Path(
    "assets/base_asset_store/gwas/ME_CFS/DecodeME/analysis"
    "/decode_mechr20_47000000_48200000_palindromes_keep_susie_finemap"
)
LD_FILE = SUSIE_DIR / "filtered_ld.npy"
GWAS_FILE = SUSIE_DIR / "filtered_gwas.parquet"
POS_COL = "POS"
SAMPLE_SECONDS = 0.05


def _rss_anon_kib() -> int:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("RssAnon:"):
            return int(line.split()[1])
    raise RuntimeError("RssAnon missing")


class _PeakSampler(threading.Thread):
    def __init__(self) -> None:
        super().__init__(daemon=True)
        self.peak_kib = 0
        self.stopped = threading.Event()

    def run(self) -> None:
        while not self.stopped.is_set():
            self.peak_kib = max(self.peak_kib, _rss_anon_kib())
            time.sleep(SAMPLE_SECONDS)


def _gib(kib: float) -> float:
    return kib / 2**20


class stage:
    """Context manager: report peak RSS and delta over the enclosed stage."""

    def __init__(self, name: str) -> None:
        self.name = name

    def __enter__(self) -> "stage":
        gc.collect()
        self.baseline = _rss_anon_kib()
        self.sampler = _PeakSampler()
        self.sampler.peak_kib = self.baseline
        print(f"\n--- STAGE {self.name}: baseline RSS {_gib(self.baseline):.2f} GiB", flush=True)
        self.sampler.start()
        return self

    def __exit__(self, *exc) -> None:
        self.sampler.stopped.set()
        self.sampler.join()
        peak = self.sampler.peak_kib
        print(
            f"--- STAGE {self.name}: peak RSS {_gib(peak):.2f} GiB "
            f"(delta {_gib(peak - self.baseline):.2f} GiB)",
            flush=True,
        )


def main() -> None:
    ld = np.load(LD_FILE)
    pos = pl.read_parquet(GWAS_FILE)[POS_COL].to_numpy()
    n_unique = int(np.unique(pos).size)
    print(
        f"ld shape={ld.shape} dtype={ld.dtype} nbytes={ld.nbytes / 2**20:.0f} MiB; "
        f"pos n={pos.size} unique={n_unique}",
        flush=True,
    )

    with stage("A square (ld**2)"):
        to_plot = ld**2

    with stage("B xarray groupby-mean (bin_options=None)"):
        ar, edges = get_array_and_edges_for_ld_heatmap(
            ld_abs=to_plot, pos=pos, bin_options=None
        )
        print(f"    binned array shape={ar.shape} edges={edges.size}", flush=True)

    with stage("C pcolormesh render + savefig"):
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.pcolormesh(edges, edges, ar, shading="auto", vmin=0, vmax=1, cmap="plasma")
        fig.savefig("/tmp/_stackplot_probe.png", dpi=100)
        plt.close(fig)

    print("\n=== profile done ===", flush=True)


if __name__ == "__main__":
    main()
