from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.analysis.bellenguez_standard_analysis import (
    BELLENGUEZ_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

BELLENGUEZ_PARQUET = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=BELLENGUEZ_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task,
    asset_id="bellenguez_keep_version_parquet_table_from_sumstats",
    sub_dir="processed",
    # pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
)


BELLENGUEZ_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "bellenguez_ppp_rg_cis_excluded",
    trait_task=BELLENGUEZ_PARQUET,
    config=PppRgConfig(variant_set="cis_excluded", trait_total_sample_size=487511),
)
