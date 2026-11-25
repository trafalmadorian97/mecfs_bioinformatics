from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_minimal_processing import (
    DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_lead_variants_task import (
    GwasLabLeadVariantsTask,
)

DECODE_ME_GWAS_1_LEAD_VARIANTS = GwasLabLeadVariantsTask(
    sumstats_task=DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
    short_id="decode_me_gwas_1_lead_variants",
)
