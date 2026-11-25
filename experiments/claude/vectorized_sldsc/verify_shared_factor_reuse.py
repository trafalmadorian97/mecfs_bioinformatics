"""
Run the vectorized shared_controls Task's execute() for every annotation set on DecodeME and
compare its output with the stored result, to check that factoring A' W A once per control
group (instead of solving it afresh for each chunk) leaves the estimates unchanged.

Usage:
    pixi r python experiments/claude/vectorized_sldsc/verify_shared_factor_reuse.py 2>&1 | tee experiments/claude/vectorized_sldsc/verify_shared_factor_reuse.log
"""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import polars as pl

from experiments.claude.vectorized_sldsc.vectorized_sldsc_experiment import _SHARED
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_COL,
    SLDSC_CTS_COEF_SE_COL,
    SLDSC_CTS_NAME_COL,
)


def main() -> None:
    for set_name, group in _SHARED.partitioned_tasks.items():
        task = group.cell_analysis_task
        assets = DEFAULT_RUNNER.run(list(task.deps))
        stored = pl.read_csv(DEFAULT_RUNNER.meta_to_path(task.meta)).sort(
            SLDSC_CTS_NAME_COL
        )
        with TemporaryDirectory() as tmp:
            result = task.execute(
                scratch_dir=Path(tmp), fetch=assets.__getitem__, wf=make_wf()
            )
            assert isinstance(result, FileAsset)
            fresh = pl.read_csv(result.path).sort(SLDSC_CTS_NAME_COL)
        assert fresh[SLDSC_CTS_NAME_COL].equals(stored[SLDSC_CTS_NAME_COL])
        rel = {
            col: np.max(
                np.abs(fresh[col].to_numpy() - stored[col].to_numpy())
                / np.abs(stored[col].to_numpy())
            )
            for col in (SLDSC_CTS_COEF_COL, SLDSC_CTS_COEF_SE_COL)
        }
        print(
            f"{set_name:30s} identical={fresh.equals(stored)} "
            f"max rel dcoef={rel[SLDSC_CTS_COEF_COL]:.2e} "
            f"max rel dse={rel[SLDSC_CTS_COEF_SE_COL]:.2e}",
            flush=True,
        )


if __name__ == "__main__":
    main()
