from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import \
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import GwasLabSumstatsToTableTask
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

JOHNSTON_ET_AL_PARQUET = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task,
    asset_id="johnston_et_al_keep_version_parquet_table_from_sumstats",
    sub_dir="processed",
    # pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
)


JOHNSTON_ET_AL_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "johnston_et_al_pain_ppp_rg_cis_excluded",
    trait_task=JOHNSTON_ET_AL_PARQUET,
    config=PppRgConfig(variant_set="cis_excluded",trait_total_sample_size=387649),
)
