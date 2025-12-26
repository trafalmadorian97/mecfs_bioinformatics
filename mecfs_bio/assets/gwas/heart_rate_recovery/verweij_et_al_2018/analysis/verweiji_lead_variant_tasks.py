from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.processed.verweiji_magma_task_generator import (
    VERWEIJI_ET_AL_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_lead_variants_task import (
    GwasLabLeadVariantsTask,
)

VERWEIJI_LEAD_VARIANTS = GwasLabLeadVariantsTask(
    sumstats_task=VERWEIJI_ET_AL_COMBINED_MAGMA_TASKS.sumstats_task,
    short_id="verweiji_lead_variants",
)
