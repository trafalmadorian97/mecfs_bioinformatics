from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_sumstats import (
    LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)

LIU_ET_AL_2023_IBD_LIFTOVER_TO_37_PARQUET = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
        asset_id="liu_et_al_2023_ibd_eur_liftover_to_37_parquet",
        sub_dir="processed",
    )
)
