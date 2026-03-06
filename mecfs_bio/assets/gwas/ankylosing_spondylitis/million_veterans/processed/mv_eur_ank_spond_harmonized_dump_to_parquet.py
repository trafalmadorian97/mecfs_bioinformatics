from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_sumstats_harmonized import (
    MILLION_VETERANS_ANK_SPOND_SUMSTATS_HARMONIZED,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)

MV_EUR_ANK_SPOND_HARMONIZED_DUMP_TO_PARQUET = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=MILLION_VETERANS_ANK_SPOND_SUMSTATS_HARMONIZED,
        asset_id="million_veterans_eur_ank_spond_harmonized_dump_to_parquet",
        sub_dir="processed",
    )
)
