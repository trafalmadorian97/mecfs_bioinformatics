"""
Run the vectorized exact Task's execute() for one annotation set on DecodeME and check its
output equals the stored result; run under /usr/bin/time -v to read peak RSS.

Usage:
    /usr/bin/time -v pixi r python experiments/claude/vectorized_sldsc/verify_drop_references.py <set_name>
"""

import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import polars as pl

from experiments.claude.vectorized_sldsc.vectorized_sldsc_experiment import _EXACT
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.wf.base_wf import make_wf


def main(set_name: str) -> None:
    task = _EXACT.partitioned_tasks[set_name].cell_analysis_task
    assets = DEFAULT_RUNNER.run(list(task.deps))
    stored = pl.read_csv(DEFAULT_RUNNER.meta_to_path(task.meta))
    with TemporaryDirectory() as tmp:
        result = task.execute(
            scratch_dir=Path(tmp), fetch=assets.__getitem__, wf=make_wf()
        )
        assert isinstance(result, FileAsset)
        fresh = pl.read_csv(result.path)
    print(f"{set_name}: identical to stored result: {fresh.equals(stored)}")


if __name__ == "__main__":
    main(sys.argv[1])
