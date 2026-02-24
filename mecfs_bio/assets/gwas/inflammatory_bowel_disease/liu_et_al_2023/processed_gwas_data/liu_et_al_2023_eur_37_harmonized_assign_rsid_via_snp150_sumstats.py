from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar import (
    LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)

LIU_ET_AL_2023_RSIDS_FROM_ANNOVAR_SUMSATS = GWASLabCreateSumstatsTask(
    df_source_task=LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
    asset_id=AssetId("liu_et_al_2023_rsids_from_annovar_sumstats"),
    basic_check=True,
    genome_build="19",
    fmt="gwaslab",
)
