"""
Task to apply S-LDSC to the aegisdottir Syncope data
"""
from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_sumstats_liftover_hapmap3_dedup import \
    AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS_HAPMAP3_DEDUP
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

SYNCOPE_S_LDSC_TASKS = standard_sldsc_task_generator(
    sumstats_task=AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS_HAPMAP3_DEDUP,
    base_name="aegisdottir_syncope_sldcs",
    pre_pipe=CompositePipe(
        [
            ComputeBetaPipe(),
            SetColToConstantPipe(
                col_name=GWASLAB_SAMPLE_SIZE_COLUMN,
                constant=946_861,
            ),
        ]
    ),
)
