"""
Task to dump summary stats to a parquet file
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_liftover_to_37 import (
    DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)

DECODE_ME_GWAS_1_LIFTOVER_TO_37_PARQUET = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
        asset_id="decode_me_gwas_1_liftover_to_37_parquet_file",
        sub_dir="processed",
    )
)
