"""
Apply MR using cis-pQTL to the IBD GWAS
"""

from mecfs_bio.asset_generator.mr_with_cis_pqtl_asset_generator import (
    BinaryOutcomeConfig,
    mr_cis_pqtl_asset_generator,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar import (
    LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
)
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

LIU_ET_AL_CIS_PQTL_MR = mr_cis_pqtl_asset_generator(
    gwas_dataframe_with_rsids=LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
    base_name="liu_et_al_2023_ibd",
    outcome_config=BinaryOutcomeConfig(
        n_case=24560,
        n_control=34_915,
    ),
    pre_pipe=RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
)
