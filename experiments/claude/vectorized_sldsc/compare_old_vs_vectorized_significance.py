"""
Compare the significant cell/tissue types of the vectorized S-LDSC runs used by the
documentation figures against those of the original (gwaslab) S-LDSC runs they replace.

For every vectorized S-LDSC task group reachable from ALL_FIGURE_TASKS, the original task group
is reconstructed with standard_sldsc_task_generator (same sumstats, same base name). Asset paths
depend only on metadata, so the stored original results are read without rebuilding anything.
Groups whose original results were never built are reported and skipped.

For each annotation set this prints the significant sets under both runs and any differences,
the largest coefficient and p-value discrepancies, and the margin by which the closest
cell type sits from the significance threshold under each run (corrected p / alpha, so a ratio
near 1 marks a fragile call).

Usage:
    pixi r python experiments/claude/vectorized_sldsc/compare_old_vs_vectorized_significance.py 2>&1 | tee experiments/claude/vectorized_sldsc/compare_old_vs_vectorized_significance.log
"""

import numpy as np
import polars as pl
from statsmodels.stats.multitest import multipletests

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.multiple_testing_table_task import (
    REJECT_NULL_LABEL,
    MultipleTestingTableTask,
)
from mecfs_bio.build_system.task.sldsc.vectorized_cell_analysis_by_sldsc_task import (
    VectorizedCellAnalysisBySLDSCTask,
)
from mecfs_bio.build_system.task_generator.sldsc_task_generator import VEC_ID_SUFFIX
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_COL,
    SLDSC_CTS_NAME_COL,
    SLDSC_CTS_P_COL,
)
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS

_CELL_SUFFIX = "_cell_analysis_by_s_ldsc" + VEC_ID_SUFFIX


def _reachable(tasks: list[Task]) -> dict[AssetId, Task]:
    seen: dict[AssetId, Task] = {}
    stack = list(tasks)
    while stack:
        task = stack.pop()
        if task.asset_id not in seen:
            seen[task.asset_id] = task
            stack.extend(task.deps)
    return seen


def _significant(task: Task) -> set[str]:
    frame = pl.read_csv(DEFAULT_RUNNER.meta_to_path(task.meta))
    if REJECT_NULL_LABEL in frame.columns:
        frame = frame.filter(pl.col(REJECT_NULL_LABEL))
    return set(frame[SLDSC_CTS_NAME_COL].to_list())


def _threshold_margin(cells: pl.DataFrame, alpha: float, method: str) -> float:
    """corrected p / alpha of the cell type closest to the threshold, on a log scale."""
    _, corrected, _, _ = multipletests(
        cells[SLDSC_CTS_P_COL].to_numpy(), alpha=alpha, method=method
    )
    return float(np.exp(np.min(np.abs(np.log(corrected / alpha)))))


def main() -> None:
    reachable = _reachable(ALL_FIGURE_TASKS)
    new_mt_tasks = sorted(
        (
            t
            for t in reachable.values()
            if isinstance(t, MultipleTestingTableTask)
            and str(t.asset_id).endswith(VEC_ID_SUFFIX)
        ),
        key=lambda t: str(t.asset_id),
    )
    n_compared = n_changed = 0
    for new_mt in new_mt_tasks:
        (new_cell,) = new_mt.deps
        assert isinstance(new_cell, VectorizedCellAnalysisBySLDSCTask)
        prefix = str(new_cell.asset_id).removesuffix(_CELL_SUFFIX)
        old_generator = standard_sldsc_task_generator(
            sumstats_task=new_cell.source_sumstats_task, base_name="__unused__"
        )
        entry = next(
            e for e in old_generator.partitioned_tasks if prefix.endswith("_" + e)
        )
        base_name = prefix.removesuffix("_" + entry)
        old_group = standard_sldsc_task_generator(
            sumstats_task=new_cell.source_sumstats_task, base_name=base_name
        ).partitioned_tasks[entry]
        old_mt, old_cell = old_group.multiple_testing_task, old_group.cell_analysis_task
        label = f"{base_name} / {entry}"
        if not (
            DEFAULT_RUNNER.meta_to_path(old_mt.meta).exists()
            and DEFAULT_RUNNER.meta_to_path(old_cell.meta).exists()
        ):
            print(f"{label}: original results not built, skipped", flush=True)
            continue
        assert isinstance(old_mt, MultipleTestingTableTask)
        n_compared += 1

        old_sig, new_sig = _significant(old_mt), _significant(new_mt)
        old_cells = pl.read_csv(DEFAULT_RUNNER.meta_to_path(old_cell.meta))
        new_cells = pl.read_csv(DEFAULT_RUNNER.meta_to_path(new_cell.meta))
        joined = old_cells.join(
            new_cells, on=SLDSC_CTS_NAME_COL, how="full", suffix="_new", coalesce=True
        )
        assert joined.height == old_cells.height == new_cells.height, (
            f"{label}: cell type names differ between runs"
        )
        coef_rel = (
            (pl.col(SLDSC_CTS_COEF_COL + "_new") - pl.col(SLDSC_CTS_COEF_COL)).abs()
            / pl.col(SLDSC_CTS_COEF_COL).abs()
        )
        p_rel = (
            (pl.col(SLDSC_CTS_P_COL + "_new") - pl.col(SLDSC_CTS_P_COL)).abs()
            / pl.col(SLDSC_CTS_P_COL)
        )
        diffs = joined.select(coef_rel.max().alias("coef"), p_rel.max().alias("p"))
        margins = (
            _threshold_margin(old_cells, alpha=old_mt.alpha, method=old_mt.method),
            _threshold_margin(new_cells, alpha=new_mt.alpha, method=new_mt.method),
        )
        same = old_sig == new_sig
        n_changed += not same
        print(
            f"{label}: {'SAME' if same else 'CHANGED'} "
            f"({len(old_sig)} old / {len(new_sig)} new significant of {old_cells.height}); "
            f"max rel diff coef={diffs['coef'][0]:.1e} p={diffs['p'][0]:.1e}; "
            f"closest-call corrected p/alpha old={margins[0]:.3f} new={margins[1]:.3f}",
            flush=True,
        )
        if not same:
            print(f"    lost: {sorted(old_sig - new_sig)}", flush=True)
            print(f"    gained: {sorted(new_sig - old_sig)}", flush=True)
    print(
        f"compared {n_compared} annotation-set runs; {n_changed} changed", flush=True
    )


if __name__ == "__main__":
    main()
