from mecfs_bio.asset_generator.mr_with_cis_pqtl_asset_generator import (
    BinaryOutcomeConfig,
    mr_cis_pqtl_asset_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

DECODE_ME_BASIC_CIS_PQTL_MR = mr_cis_pqtl_asset_generator(
    gwas_dataframe_with_rsids=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    base_name="decode_me_gwas_1_mr",
    outcome_config=BinaryOutcomeConfig(n_case=15_579, n_control=259_909),
    pre_pipe=RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
    significant_alpha=0.05,
)
