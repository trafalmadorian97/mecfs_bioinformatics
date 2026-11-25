"""
Vectorized S-LDSC validation on DecodeME.

Modes (one per process so peak RSS is per configuration):
    profile-old   cProfile the existing CellAnalysisByLDSCTask on GTEx_brain (13 models), to see
                  whether time goes to file parsing or to Hsq (ldsc's "h-squared" heritability
                  regression class). Also saves that fresh gwaslab result as a second oracle.
    run-exact     build prepared references (first run only) and vectorized exact results.
    run-shared    vectorized shared_controls results.
    time-exact    per-trait cost once prepared references exist: execute each vectorized exact
                  Task directly (no store write) against the cached prepared references.
    time-shared   the same for shared_controls.
    compare       compare old (cached gwaslab) vs exact vs shared for all six sets, plus
                  exact vs the fresh gwaslab GTEx_brain result from profile-old.

The cached old results are read straight from the asset store (located via meta_to_path), so
that a stale trace can never trigger an hours-long rebuild of the old Task.

Usage:
    pixi r python experiments/claude/vectorized_sldsc/vectorized_sldsc_experiment.py <mode> 2>&1 | tee experiments/claude/vectorized_sldsc/<mode>.log
"""

import cProfile
import pstats
import resource
import shutil
import sys
import time
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Mapping, Protocol

import numpy as np
import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.concrete_vectorized_sldsc_generator import (
    standard_vectorized_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import (
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
    DECODE_ME_S_LDSC,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePIfNeededPipe
from mecfs_bio.build_system.task_generator.sldsc_task_generator import (
    CellAnalysisTaskGroup,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_COL,
    SLDSC_CTS_COEF_SE_COL,
    SLDSC_CTS_NAME_COL,
    SLDSC_CTS_P_COL,
)

_FDR_ALPHA = 0.01
_PROFILED_SET = "gtex_brain"
_FRESH_OLD_RESULT = Path(__file__).parent / "fresh_gwaslab_gtex_brain.csv"
_EXACT = standard_vectorized_sldsc_task_generator(
    sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
    base_name="decode_me_gwas_1",
    pre_pipe=ComputePIfNeededPipe(),
)
_SHARED = standard_vectorized_sldsc_task_generator(
    sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
    base_name="decode_me_gwas_1_shared_controls",
    pre_pipe=ComputePIfNeededPipe(),
    weighting="shared_controls",
)


class _CellAnalysisGenerator(Protocol):
    @property
    def partitioned_tasks(self) -> Mapping[str, CellAnalysisTaskGroup]: ...


def _peak_rss_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2


def _cell_tasks(generator: _CellAnalysisGenerator) -> dict[str, Task]:
    return {
        name: group.cell_analysis_task
        for name, group in generator.partitioned_tasks.items()
    }


def profile_old() -> None:
    task = _cell_tasks(DECODE_ME_S_LDSC)[_PROFILED_SET]
    assets = DEFAULT_RUNNER.run(list(task.deps))
    with TemporaryDirectory() as tmp:
        profiler = cProfile.Profile()
        start = time.perf_counter()
        result = profiler.runcall(
            task.execute, scratch_dir=Path(tmp), fetch=assets.__getitem__, wf=make_wf()
        )
        print(f"old {_PROFILED_SET} wall time: {time.perf_counter() - start:.1f}s")
        assert isinstance(result, FileAsset)
        shutil.copyfile(result.path, _FRESH_OLD_RESULT)
    pstats.Stats(profiler).sort_stats("cumulative").print_stats(40)
    pstats.Stats(profiler).sort_stats("tottime").print_stats(20)
    print(f"peak RSS {_peak_rss_gb():.2f} GB")


def run(generator: _CellAnalysisGenerator, label: str) -> None:
    for name, task in _cell_tasks(generator).items():
        start = time.perf_counter()
        DEFAULT_RUNNER.run([task])
        print(
            f"{label} {name}: {time.perf_counter() - start:.1f}s "
            "(includes prepared-reference build on first run)"
        )
    print(f"{label} peak RSS {_peak_rss_gb():.2f} GB")


def time_direct(generator: _CellAnalysisGenerator, label: str) -> None:
    tasks = _cell_tasks(generator)
    assets = DEFAULT_RUNNER.run([dep for task in tasks.values() for dep in task.deps])
    for name, task in tasks.items():
        with TemporaryDirectory() as tmp:
            start = time.perf_counter()
            task.execute(scratch_dir=Path(tmp), fetch=assets.__getitem__, wf=make_wf())
            print(f"{label} {name}: {time.perf_counter() - start:.1f}s")
    print(f"{label} peak RSS {_peak_rss_gb():.2f} GB")


def _bh_reject(p: np.ndarray, alpha: float) -> np.ndarray:
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    passed = np.zeros(len(p), dtype=bool)
    below = np.nonzero(ranked <= alpha)[0]
    if below.size:
        passed[order[: below.max() + 1]] = True
    return passed


def _read_stored(task: Task) -> pl.DataFrame:
    path = DEFAULT_RUNNER.meta_to_path(task.meta)
    assert path.is_file(), f"{task.asset_id} is not in the asset store at {path}"
    return pl.read_csv(path).sort(SLDSC_CTS_NAME_COL)


def _report(name: str, label: str, old: pl.DataFrame, new: pl.DataFrame) -> None:
    assert new[SLDSC_CTS_NAME_COL].to_list() == old[SLDSC_CTS_NAME_COL].to_list()
    old_se = old[SLDSC_CTS_COEF_SE_COL].to_numpy()
    z_old = old[SLDSC_CTS_COEF_COL].to_numpy() / old_se
    z_new = (new[SLDSC_CTS_COEF_COL] / new[SLDSC_CTS_COEF_SE_COL]).to_numpy()
    dz = np.abs(z_new - z_old)
    dcoef = (
        np.abs(new[SLDSC_CTS_COEF_COL].to_numpy() - old[SLDSC_CTS_COEF_COL].to_numpy())
        / old_se
    )
    fdr_changes = int(
        np.sum(
            _bh_reject(new[SLDSC_CTS_P_COL].to_numpy(), _FDR_ALPHA)
            != _bh_reject(old[SLDSC_CTS_P_COL].to_numpy(), _FDR_ALPHA)
        )
    )
    print(
        f"{name:30s} {label:14s} n={old.height:4d} max|dz|={dz.max():.2e} "
        f"median|dz|={np.median(dz):.2e} max|dcoef|/se={dcoef.max():.2e} "
        f"FDR changes={fdr_changes}"
    )


def compare() -> None:
    old_tasks = _cell_tasks(DECODE_ME_S_LDSC)
    for name, old_task in old_tasks.items():
        old = _read_stored(old_task)
        for label, generator in (("exact", _EXACT), ("shared", _SHARED)):
            _report(
                name=name,
                label=label,
                old=old,
                new=_read_stored(_cell_tasks(generator)[name]),
            )
    if _FRESH_OLD_RESULT.is_file():
        fresh = pl.read_csv(_FRESH_OLD_RESULT).sort(SLDSC_CTS_NAME_COL)
        _report(
            name=_PROFILED_SET,
            label="exact-vs-fresh",
            old=fresh,
            new=_read_stored(_cell_tasks(_EXACT)[_PROFILED_SET]),
        )
        _report(
            name=_PROFILED_SET,
            label="cached-vs-fresh",
            old=fresh,
            new=_read_stored(old_tasks[_PROFILED_SET]),
        )


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "profile-old":
        profile_old()
    elif mode == "run-exact":
        run(generator=_EXACT, label="exact")
    elif mode == "run-shared":
        run(generator=_SHARED, label="shared")
    elif mode == "time-exact":
        time_direct(generator=_EXACT, label="exact")
    elif mode == "time-shared":
        time_direct(generator=_SHARED, label="shared")
    elif mode == "compare":
        compare()
    else:
        raise ValueError(f"unknown mode {mode}")
