"""
Million Veteran Program EUR ankylosing spondylitis sumstats, genome-reference harmonized to build 37, loaded as
gwaslab sumstats for gwaslab-based analyses.
"""

from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_genome_reference_harmonized import (
    MV_EUR_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)

MILLION_VETERANS_ANK_SPOND_SUMSTATS_HARMONIZED = GWASLabCreateSumstatsTask(
    df_source_task=MV_EUR_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
    target_asset_id=AssetId(
        "million_veterans_ank_spond_genome_reference_harmonized_sumstats"
    ),
    basic_check=True,
    genome_build="19",
    fmt="gwaslab",
)
