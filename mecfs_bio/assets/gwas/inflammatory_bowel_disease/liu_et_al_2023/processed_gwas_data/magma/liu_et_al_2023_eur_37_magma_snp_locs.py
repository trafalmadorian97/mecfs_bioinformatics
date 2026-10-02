from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_with_rsid import (
    LIU_ET_AL_2023_IBD_LIFTOVER_TO_37_WITH_RSID,
)
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask

LIU_ET_AL_IBD_EUR_BUILD_37_MAGMA_SNP_LOCS = (
    MagmaSNPFileTask.create_for_magma_snp_pos_file(
        gwas_parquet_with_rsids_task=LIU_ET_AL_2023_IBD_LIFTOVER_TO_37_WITH_RSID,
        asset_id="liu_et_al_2023_ibd_eur_build_37_magma_snp_locs",
    )
)
