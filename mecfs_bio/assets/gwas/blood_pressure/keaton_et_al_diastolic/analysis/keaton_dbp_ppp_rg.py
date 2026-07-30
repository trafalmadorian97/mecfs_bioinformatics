from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import (
    KEATON_DBP_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

KEATON_DBP_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "keaton_dbp_ppp_rg_cis_excluded",
    trait_task=KEATON_DBP_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    config=PppRgConfig(variant_set="cis_excluded"),
)
