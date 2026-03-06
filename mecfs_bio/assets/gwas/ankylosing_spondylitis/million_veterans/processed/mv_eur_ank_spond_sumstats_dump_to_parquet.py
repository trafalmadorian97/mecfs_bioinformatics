from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_sumstats import (
    MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_SUMSTATS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)

MILLION_VETERANS_ANK_SPOND_SUMSTATS_37_DUMP_TO_PARQUET = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_SUMSTATS,
        asset_id="million_veterans_ank_spond_sumstats_37_dump_to_parquet",
        sub_dir="processed",
    )
)
