from mecfs_bio.assets.gwas.ankylosing_spondylitis.multi_study.fixed_effects_meta_analysis.processed.ank_spond_meta_analysis import (
    ANK_SPOND_FIXED_EFFECTS_META_ANALYSIS,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)

ANK_SPOND_META_ANALYSIS_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=ANK_SPOND_FIXED_EFFECTS_META_ANALYSIS,
    target_asset_id=AssetId("ank_spond_meta_analysis_sumstats"),
    basic_check=True,
    genome_build="19",
    fmt="gwaslab",
)
