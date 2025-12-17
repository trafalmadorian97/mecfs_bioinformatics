"""
Pipe to apply S-LSDC to the Million Veterans GWAS of red blood cell volume
"""

from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    standard_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.red_blood_cell_volume.million_veterans.analysis.rbc_volume_magma_task_generator import (
    MILLION_VETERANS_EUR_RBC_VOLUME_MAGMA_TASKS,
)

MILLION_VETERANS_EUR_RBC_VOL_S_LDSC_TASKS = standard_sldsc_task_generator(
    sumstats_task=MILLION_VETERANS_EUR_RBC_VOLUME_MAGMA_TASKS.sumstats_task,
    base_name="million_veterans_eur_rbc_vol",
    # pre_pipe=SetColToConstantPipe(
    #     col_name=GWASLAB_SAMPLE_SIZE_COLUMN,
    #     constant=407294 , # source: https://www.ebi.ac.uk/gwas/studies/GCST90475466
    # ),
)
