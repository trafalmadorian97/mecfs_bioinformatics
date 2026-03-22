from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID = GWASLabCreateSumstatsTask(
    df_source_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task,
    asset_id=AssetId("decode_me_sumstats_37_with_annovar_rsids"),
    basic_check=True,  # already checked
    genome_build="19",
    pre_pipe=RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
    fmt="gwaslab",
)
