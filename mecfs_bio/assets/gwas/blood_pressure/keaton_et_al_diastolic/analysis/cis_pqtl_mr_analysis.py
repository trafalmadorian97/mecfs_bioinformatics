from mecfs_bio.asset_generator.mr_with_cis_pqtl_asset_generator import mr_cis_pqtl_asset_generator
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import \
    KEATON_DBP_STANDARD_ANALYSIS

KEATON_DBP_CIS_PQTL_MR_ANALYSIS=mr_cis_pqtl_asset_generator(
    gwas_dataframe_with_rsids=KEATON_DBP_STANDARD_ANALYSIS.magma_tasks.parquet_file_task
)