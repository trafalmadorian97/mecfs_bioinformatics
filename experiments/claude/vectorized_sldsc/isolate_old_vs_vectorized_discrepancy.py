"""
Isolate why a vectorized S-LDSC result differs from the stored original (gwaslab) result for one
standard-analysis trait and annotation set.

Recomputes, in a temporary directory (nothing is written to the asset store), the original
gwaslab Task and the vectorized Task, each with and without the duplicate-rsID drop that the
vectorized standard generator adds, and compares every variant with the stored original and
stored vectorized results (max relative difference in coefficient and p-value).

Usage:
    pixi r python experiments/claude/vectorized_sldsc/isolate_old_vs_vectorized_discrepancy.py <base_name> <entry> 2>&1 | tee experiments/claude/vectorized_sldsc/isolate_<base_name>_<entry>.log
"""

import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import attrs
import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.sldsc.vectorized_cell_analysis_by_sldsc_task import (
    VectorizedCellAnalysisBySLDSCTask,
)
from mecfs_bio.build_system.task_generator.sldsc_task_generator import VEC_ID_SUFFIX
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_COL,
    SLDSC_CTS_NAME_COL,
    SLDSC_CTS_P_COL,
)
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS


def _find_vectorized_cell_task(cell_id: str) -> VectorizedCellAnalysisBySLDSCTask:
    stack: list[Task] = list(ALL_FIGURE_TASKS)
    seen: set[str] = set()
    while stack:
        task = stack.pop()
        if str(task.asset_id) in seen:
            continue
        seen.add(str(task.asset_id))
        if str(task.asset_id) == cell_id:
            assert isinstance(task, VectorizedCellAnalysisBySLDSCTask)
            return task
        stack.extend(task.deps)
    raise ValueError(f"{cell_id} is not reachable from ALL_FIGURE_TASKS")


def _execute(task: Task) -> pl.DataFrame:
    assets = DEFAULT_RUNNER.run(list(task.deps))
    with TemporaryDirectory() as tmp:
        result = task.execute(scratch_dir=Path(tmp), fetch=assets.__getitem__, wf=make_wf())
        assert isinstance(result, FileAsset)
        return pl.read_csv(result.path)


def _stored(task: Task) -> pl.DataFrame:
    return pl.read_csv(DEFAULT_RUNNER.meta_to_path(task.meta))


def _max_rel_diff(a: pl.DataFrame, b: pl.DataFrame) -> str:
    joined = a.join(b, on=SLDSC_CTS_NAME_COL, how="inner", suffix="_b")
    assert joined.height == a.height == b.height
    out = []
    for col in (SLDSC_CTS_COEF_COL, SLDSC_CTS_P_COL):
        rel = ((joined[col + "_b"] - joined[col]).abs() / joined[col].abs()).max()
        out.append(f"{col}={rel:.1e}")
    return " ".join(out)


def main() -> None:
    base_name, entry = sys.argv[1], sys.argv[2]
    vectorized = _find_vectorized_cell_task(
        f"{base_name}_{entry}_cell_analysis_by_s_ldsc{VEC_ID_SUFFIX}"
    )
    with_dedup = vectorized.pre_pipe
    assert isinstance(with_dedup, CompositePipe)
    without_dedup = with_dedup.pipes[0]
    original = standard_sldsc_task_generator(
        sumstats_task=vectorized.source_sumstats_task,
        base_name=base_name,
        pre_pipe=without_dedup,
    ).partitioned_tasks[entry].cell_analysis_task

    results = {
        "stored original": _stored(original),
        "stored vectorized": _stored(vectorized),
    }
    runs = {
        "fresh original": original,
        "fresh original + dedup": attrs.evolve(original, pre_pipe=with_dedup),
        "fresh vectorized - dedup": attrs.evolve(vectorized, pre_pipe=without_dedup),
        "fresh vectorized": vectorized,
    }
    for label, task in runs.items():
        print(f"running {label}...", flush=True)
        results[label] = _execute(task)
    labels = list(results)
    for i, a in enumerate(labels):
        for b in labels[i + 1 :]:
            print(f"{a:26s} vs {b:26s}: {_max_rel_diff(results[a], results[b])}", flush=True)


if __name__ == "__main__":
    main()
