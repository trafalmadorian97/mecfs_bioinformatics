from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.processed.finngen_ank_spond_harmonized_dump_to_parquet import \
    FINNGEN_ANK_SPOND_HARMONIZED_DUMP_TO_PARQUET
from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_harmonized_dump_to_parquet import \
    MV_EUR_ANK_SPOND_HARMONIZED_DUMP_TO_PARQUET
from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_ank_spond_harmonized_dump_to_parquet import \
    UKBB_ANK_SPOND_HARMONIZED_DUMP_TO_PARQUET
from mecfs_bio.build_system.task.fixed_effect_meta_analysis_task import FixedEffectsMetaAnalysisTask, GwasSource, \
    CaseControlSampleInfo
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, \
    GWASLAB_NON_EFFECT_ALLELE_COL, GWASLAB_SE_COL

ANK_SPOND_FIXED_EFFECTS_META_ANALYSIS=FixedEffectsMetaAnalysisTask.create(
    asset_id="ank_spond_fixed_effects_meta_analysis",
    meta_analysis_name="fixed_effects_meta_analysis",
    sources=[
        GwasSource(
            task=FINNGEN_ANK_SPOND_HARMONIZED_DUMP_TO_PARQUET,
            sample_info=CaseControlSampleInfo(
                cases=1462,
                controls=164682
            ),#https://opengwas.io/datasets/finn-b-M13_ANKYLOSPON#
            pipe=UniquePipe( # Most the duplications in this gwas appear to be due to anomalies caused by conversion of multi-allelic to bi-allelic variants
                by=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL
                ],
                keep="none",
                order_by=[GWASLAB_POS_COL]
            )
        ),
        GwasSource(
            task=MV_EUR_ANK_SPOND_HARMONIZED_DUMP_TO_PARQUET,
            sample_info=CaseControlSampleInfo(
                cases=1637,
                controls=450630, # source:https://www.ebi.ac.uk/gwas/studies/GCST90476232
            ),

            pipe=CompositePipe([

                ComputeBetaPipe(),
                ComputeSEPipe(),

                UniquePipe(
                # duplicated rows in the million veterans GWAS are likely the consequence of liftover effects
                by=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL
                ],
                keep="first",
                order_by=[GWASLAB_SE_COL]

            ),


            ],

            )

        ),
        GwasSource(
            task=UKBB_ANK_SPOND_HARMONIZED_DUMP_TO_PARQUET,
            sample_info=CaseControlSampleInfo(
                    cases=2076,
                controls=456364, # https://www.ebi.ac.uk/gwas/studies/GCST90474065
            ),

            pipe=CompositePipe([UniquePipe(
                # duplicated rows in the million veterans GWAS are likely the consequence of liftover effects
                by=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL
                ],
                keep="first",
                order_by=[GWASLAB_SE_COL]

            ),

                # ComputeBetaPipe(),
                # ComputeSEPipe(),

            ],

            )

        )
    ]
)