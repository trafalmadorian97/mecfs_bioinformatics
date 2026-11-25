"""
Real-bucket smoke run of the deCODE SomaScan database tasks (plan Task 17, step 1), one stage
per process so /usr/bin/time -v gives each stage's own peak RSS:

  pixi r /usr/bin/time -v python experiments/claude/decode_proteomics_survey/smoke_run.py <stage>

stage: variant_info | clean_variant_info | index | protein:<seq_id>
"""

import sys
import time

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_database_protein_files import (
    HAPMAP_3_DECODE_SMP_DATABASE,
)
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_index import (
    HAPMAP_3_DECODE_SMP_INDEX,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_variant_info import (
    DECODE_CLEAN_VARIANT_INFO,
    DECODE_VARIANT_INFO,
)
from mecfs_bio.build_system.task.base_task import Task


def _target(stage: str) -> Task:
    if stage == "variant_info":
        return DECODE_VARIANT_INFO
    if stage == "clean_variant_info":
        return DECODE_CLEAN_VARIANT_INFO
    if stage == "index":
        return HAPMAP_3_DECODE_SMP_INDEX
    seq_id = stage.removeprefix("protein:")
    (task,) = [
        t
        for t in HAPMAP_3_DECODE_SMP_DATABASE.protein_tasks
        if t.protein.seq_id == seq_id
    ]
    return task


def main() -> None:
    target = _target(sys.argv[1])
    start = time.perf_counter()
    assets = DEFAULT_RUNNER.run([target], incremental_save=True)
    print(
        f"stage {sys.argv[1]}: {time.perf_counter() - start:.1f}s -> {assets[target.asset_id]}"
    )


if __name__ == "__main__":
    main()
