from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.processed.finngen_ank_spond_sumstats_harmonized import (
    FINGNEN_ANK_SPOND_SUMSTATS_HARMONIZED,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)

FINNGEN_ANK_SPOND_HARMONIZED_DUMP_TO_PARQUET = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=FINGNEN_ANK_SPOND_SUMSTATS_HARMONIZED,
        asset_id="finngen_ank_spond_harmonized_dump_to_parquet",
        sub_dir="processed",
    )
)
