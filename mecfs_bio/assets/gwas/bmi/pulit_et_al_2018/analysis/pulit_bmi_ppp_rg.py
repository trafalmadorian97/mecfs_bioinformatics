"""
Genetic correlation between the GIANT + UK Biobank GWAS of body mass index (Pulit et al. 2019) and
plasma protein levels from the UK Biobank Pharma Proteomics Project, with each protein's cis window
excluded.

The BMI sumstats carry a per-variant N column, so no constant trait sample size is supplied.
"""

from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.analysis.pulit_bmi_standard_analysis import (
    PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

PULIT_ET_AL_2018_BMI_PARQUET = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task,
    asset_id="pulit_et_al_2018_bmi_keep_version_parquet_table_from_sumstats",
    sub_dir="processed",
)


PULIT_ET_AL_2018_BMI_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "pulit_et_al_2018_bmi_ppp_rg_cis_excluded",
    trait_task=PULIT_ET_AL_2018_BMI_PARQUET,
    config=PppRgConfig(variant_set="cis_excluded"),
)
