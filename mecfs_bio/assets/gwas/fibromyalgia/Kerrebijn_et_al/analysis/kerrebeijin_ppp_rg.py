from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.analysis.standard_analysis_kerrebijin_fibro import (
    KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

KERREBEIJIN_ET_AL_PARQUET = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task,
    asset_id="kerrebeijin_et_al_keep_version_parquet_table_from_sumstats",
    sub_dir="processed",
)


KERREBEIJIN_ET_AL_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "kerrebeijin_et_al_pain_ppp_rg_cis_excluded",
    trait_task=KERREBEIJIN_ET_AL_PARQUET,
    config=PppRgConfig(variant_set="cis_excluded"),
)
