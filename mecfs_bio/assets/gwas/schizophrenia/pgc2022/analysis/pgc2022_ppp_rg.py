from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import (
    SCH_PGC_2022_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

PGC2022_SCH_PARQUET = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=SCH_PGC_2022_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task,
    asset_id="pfc_2022_keep_version_parquet_table_from_sumstats",
    sub_dir="processed",
    # pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
)


PGC2022_SCH_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "pgc2022_sch_ppp_rg_cis_excluded",
    trait_task=PGC2022_SCH_PARQUET,
    config=PppRgConfig(variant_set="cis_excluded"),
)
