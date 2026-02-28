from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_sumstats_liftover import \
    AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import GwasLabSumstatsToTableTask

AEGISDOTTIR_SUMSTATS_DUMP_TO_PARQUET =GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS,
    asset_id="aegissdottir_sumstats_37_dump_to_parquet_task",
sub_dir = "processed",
)