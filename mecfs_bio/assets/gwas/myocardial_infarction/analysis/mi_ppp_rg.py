from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_standard_analysis import \
    MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import GwasLabSumstatsToTableTask
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

MV_MI_PARQUET= GwasLabSumstatsToTableTask.create_from_source_task(
            source_tsk=MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task,
            asset_id="mv_mi_keep_version_parquet_table_from_sumstats",
            sub_dir="processed",
    pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()])
        )


MV_MI_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "mv_myocardial_infarction_ppp_rg_cis_excluded",
    trait_task=MV_MI_PARQUET,
    config=PppRgConfig(variant_set="cis_excluded",trait_total_sample_size=432053),
)
