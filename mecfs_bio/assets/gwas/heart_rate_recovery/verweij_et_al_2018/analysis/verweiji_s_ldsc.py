from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.processed.verweiji_magma_task_generator import (
    VERWEIJI_ET_AL_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import (
    GWASLAB_SAMPLE_SIZE_COLUMN,
)
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe

# todo:
# Try again with a different info filtering strategy
VERWEIJI_SLDSC_TASK_GROUP = standard_sldsc_task_generator(
    sumstats_task=VERWEIJI_ET_AL_COMBINED_MAGMA_TASKS.sumstats_task,
    base_name="verweiji_et_al_2023",
    pre_pipe=SetColToConstantPipe(
        col_name=GWASLAB_SAMPLE_SIZE_COLUMN,
        constant=58_818,  # source: abstract
    ),
)
