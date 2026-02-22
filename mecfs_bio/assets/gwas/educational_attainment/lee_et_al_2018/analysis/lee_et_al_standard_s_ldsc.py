from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import (
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_SAMPLE_SIZE_COLUMN,
)

LEE_ET_AL_EDU_STANDARD_SLDSC_TASK_GROUP = standard_sldsc_task_generator(
    sumstats_task=LEE_ET_AL_2018_COMBINED_MAGMA_TASKS.sumstats_task,
    base_name="lee_et_al_2018",
    pre_pipe=SetColToConstantPipe(
        col_name=GWASLAB_SAMPLE_SIZE_COLUMN,
        constant=257841,  # source: 30038396-GCST006572-EFO_0008354.h.tsv.gz-meta.yaml
    ),
)
