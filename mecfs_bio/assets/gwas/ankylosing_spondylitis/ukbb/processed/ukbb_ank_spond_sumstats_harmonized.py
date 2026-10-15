"""
UK Biobank ankylosing spondylitis sumstats, genome-reference harmonized to build 37, loaded as
gwaslab sumstats for gwaslab-based analyses.
"""

from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_ank_spond_genome_reference_harmonized import (
    UKBB_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)

UKBB_ANK_SPOND_SUMSTATS_HARMONIZED = GWASLabCreateSumstatsTask(
    df_source_task=UKBB_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
    target_asset_id=AssetId(
        "uk_biobank_ank_spond_genome_reference_harmonized_sumstats"
    ),
    basic_check=True,
    genome_build="19",
    fmt="gwaslab",
)
