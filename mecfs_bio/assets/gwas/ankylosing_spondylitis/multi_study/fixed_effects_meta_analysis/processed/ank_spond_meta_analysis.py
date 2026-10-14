from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.processed.finngen_ank_spond_genome_reference_harmonized import (
    FINNGEN_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_genome_reference_harmonized import (
    MV_EUR_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_ank_spond_genome_reference_harmonized import (
    UKBB_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
)
from mecfs_bio.build_system.task.fixed_effect_meta_analysis_task import (
    CaseControlSampleInfo,
    FixedEffectsMetaAnalysisTask,
    GwasSource,
)

ANK_SPOND_FIXED_EFFECTS_META_ANALYSIS = FixedEffectsMetaAnalysisTask.create(
    asset_id="ank_spond_fixed_effects_meta_analysis",
    meta_analysis_name="fixed_effects_meta_analysis",
    sources=[
        GwasSource(
            task=FINNGEN_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
            sample_info=CaseControlSampleInfo(
                cases=1462, controls=164682
            ),  # https://opengwas.io/datasets/finn-b-M13_ANKYLOSPON#
        ),
        GwasSource(
            task=MV_EUR_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
            sample_info=CaseControlSampleInfo(
                cases=1637,
                controls=450630,  # source:https://www.ebi.ac.uk/gwas/studies/GCST90476232
            ),
        ),
        GwasSource(
            task=UKBB_ANK_SPOND_GENOME_REFERENCE_HARMONIZED,
            sample_info=CaseControlSampleInfo(
                cases=2076,
                controls=456364,  # https://www.ebi.ac.uk/gwas/studies/GCST90474065
            ),
        ),
    ],
)
