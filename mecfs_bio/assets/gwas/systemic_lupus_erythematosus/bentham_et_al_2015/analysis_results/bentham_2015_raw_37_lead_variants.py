from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.processed_gwas_data.bentham_2015_raw_sumstats_minimal_processing import (
    BENTHAM_2015_SUMSTATS_MINIMAL_PROCESSING_FROM_RAW_BUILD_37,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_lead_variants_task import (
    GwasLabLeadVariantsTask,
)

BENTHAM_2015_RAW_LEAD_VARIANTS = GwasLabLeadVariantsTask(
    sumstats_task=BENTHAM_2015_SUMSTATS_MINIMAL_PROCESSING_FROM_RAW_BUILD_37,
    short_id="bentham_2015_raw_lead_variants",
)
