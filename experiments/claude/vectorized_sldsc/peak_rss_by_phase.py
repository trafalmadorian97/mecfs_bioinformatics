"""
Peak RSS of each phase of VectorizedCellAnalysisBySLDSCTask.execute on DecodeME, for one
annotation set, to find where the vectorized S-LDSC memory goes.

A background thread samples VmRSS from /proc every 20 ms; each phase reports its own peak and
the RSS left when it ends. The phases replicate execute() step by step. With keep, the same
objects stay alive as in execute(); with drop, the full prepared annotation set and baseline are
released before the kernel (only the row-aligned copies are needed there).

Usage:
    pixi r python experiments/claude/vectorized_sldsc/peak_rss_by_phase.py <set_name> <weighting> <chunk_size> <keep|drop> 2>&1 | tee experiments/claude/vectorized_sldsc/peak_rss_<set_name>_<weighting>.log
"""

import gc
import sys
import threading
import time
from pathlib import Path

from experiments.claude.vectorized_sldsc.vectorized_sldsc_experiment import (
    _EXACT,
    _SHARED,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.task.sldsc.batched_sldsc_cts import batched_sldsc_cts
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    read_prepared_annotation_set,
    read_prepared_baseline,
)
from mecfs_bio.build_system.task.sldsc.sldsc_regression_inputs import (
    build_regression_inputs,
    load_sldsc_sumstats,
)
from mecfs_bio.build_system.task.sldsc.vectorized_cell_analysis_by_sldsc_task import (
    VectorizedCellAnalysisBySLDSCTask,
)
from mecfs_bio.constants.sldsc_constants import SLDSCWeighting

_SAMPLE_SECONDS = 0.02


def _rss_gb() -> float:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) / 1024**2
    raise AssertionError("no VmRSS")


class _PeakSampler:
    def __init__(self) -> None:
        self.peak = _rss_gb()
        self._stop = False
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()

    def _loop(self) -> None:
        while not self._stop:
            self.peak = max(self.peak, _rss_gb())
            time.sleep(_SAMPLE_SECONDS)

    def phase(self, name: str, start: float) -> None:
        self.peak = max(self.peak, _rss_gb())
        print(
            f"{name:28s} {time.perf_counter() - start:7.1f}s  "
            f"peak {self.peak:6.2f} GB  end {_rss_gb():6.2f} GB",
            flush=True,
        )
        self.peak = _rss_gb()


def main(
    set_name: str, weighting: SLDSCWeighting, chunk_size: int, drop_references: bool
) -> None:
    generator = _EXACT if weighting == "exact" else _SHARED
    task = generator.partitioned_tasks[set_name].cell_analysis_task
    assert isinstance(task, VectorizedCellAnalysisBySLDSCTask)
    assets = DEFAULT_RUNNER.run(list(task.deps))
    sampler = _PeakSampler()
    start = time.perf_counter()
    sampler.phase("startup (imports, runner)", start)

    start = time.perf_counter()
    sumstats = load_sldsc_sumstats(
        asset=assets[task.source_sumstats_task.asset_id], pre_pipe=task.pre_pipe
    )
    sampler.phase("load_sldsc_sumstats", start)

    start = time.perf_counter()
    baseline_asset = assets[task.baseline_task.asset_id]
    assert isinstance(baseline_asset, DirectoryAsset)
    baseline = read_prepared_baseline(baseline_asset.path)
    sampler.phase("read_prepared_baseline", start)

    start = time.perf_counter()
    annotation_asset = assets[task.annotation_set_task.asset_id]
    assert isinstance(annotation_asset, DirectoryAsset)
    annotation_set = read_prepared_annotation_set(annotation_asset.path)
    print(
        f"  annotation matrix {annotation_set.ld.matrix.shape} "
        f"{annotation_set.ld.matrix.nbytes / 1024**3:.2f} GB"
    )
    sampler.phase("read_prepared_annotation_set", start)

    start = time.perf_counter()
    aligned = build_regression_inputs(
        sumstats=sumstats,
        baseline=baseline,
        annotation_set=annotation_set,
        min_hapmap3_snps=task.min_hapmap3_snps,
        duplicate_rsids=task.duplicate_rsids,
    )
    print(
        f"  aligned annotation_ld {aligned.inputs.annotation_ld.shape} "
        f"{aligned.inputs.annotation_ld.nbytes / 1024**3:.2f} GB"
    )
    sampler.phase("build_regression_inputs", start)

    if drop_references:
        start = time.perf_counter()
        del annotation_set, baseline, sumstats
        gc.collect()
        sampler.phase("drop full references", start)

    start = time.perf_counter()
    batched_sldsc_cts(
        inputs=aligned.inputs,
        models=aligned.models,
        weighting=task.weighting,
        chunk_size=chunk_size,
    )
    sampler.phase(f"kernel (chunk_size={chunk_size})", start)


if __name__ == "__main__":
    weighting_arg = sys.argv[2]
    assert weighting_arg in ("exact", "shared_controls")
    main(
        set_name=sys.argv[1],
        weighting=weighting_arg,
        chunk_size=int(sys.argv[3]),
        drop_references={"keep": False, "drop": True}[sys.argv[4]],
    )
