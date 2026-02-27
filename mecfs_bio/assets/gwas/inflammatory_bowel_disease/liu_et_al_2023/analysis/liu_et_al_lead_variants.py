from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_with_rsid_sumstats import (
    LIU_ET_AL_SUMSTATS_WITH_RSID_FROM_SNP150,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_lead_variants_task import (
    GwasLabLeadVariantsTask,
)

LIU_ET_AL_LEAD_VARIANTS = GwasLabLeadVariantsTask(
    short_id=AssetId("liu_et_al_2023_ibd_lead_variants"),
    sumstats_task=LIU_ET_AL_SUMSTATS_WITH_RSID_FROM_SNP150,
)
