"""
Task to use cross-trait LD Score Regression to estimate genetic correlation between a several traits.
"""

from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.asthma_standard_analysis import (
    HAN_ASTHMA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import (
    KEATON_DBP_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_with_annovar_37_rsids_sumstats import (
    DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import (
    SCH_PGC_2022_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    GeneticCorrelationByCTLDSCTask,
    QuantPhenotype,
    SumstatsSource,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

CT_LDSC_INITIAL = GeneticCorrelationByCTLDSCTask.create(
    "initial_genetic_correlation_by_ct_ldsc",
    sources=[
        SumstatsSource(
            DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
            alias="DecodeME",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                constant=275488,
            ),  # true total sample size. From preprint
            sample_info=BinaryPhenotypeSampleInfo(
                sample_prevalence=0.0566,  # 15,579/(259,909+15,579).  See: https://www.medrxiv.org/content/10.1101/2025.08.06.25333109v1.full-text#:~:text=Our%20primary%20GWAS%2C%20GWAS%2D1%2C%20compared%2015%2C579%20DecodeME%20cases%20with%20259%2C909%20UKB%20controls%20(case%3Acontrol%20ratio%20of%201%3A17)%2C%20across%20all%20autosomes
                estimated_population_prevalence=0.006,  # See Samms and Ponting: https://pmc.ncbi.nlm.nih.gov/articles/PMC12120426/#:~:text=our%20estimated%20ME/CFS%20prevalence%20in%20the%20UK%20then%20rises%20from%20330%2C000%20to%20410%2C000%20(0.6%25)
            ),
        ),
        SumstatsSource(
            JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="multisite_pain",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN, constant=387649
            ),  # True total sample size. From Gwas catalog
            sample_info=QuantPhenotype(),
        ),
        SumstatsSource(
            SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="schizophrenia",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                constant=130644,  # from summary statistics file (130644=53386+77258)
            ),  # True total sample size. From Gwas catalog
            sample_info=BinaryPhenotypeSampleInfo(
                sample_prevalence=0.408,  # 53386/(53386 +77258)
                estimated_population_prevalence=0.01,  # See: https://pmc.ncbi.nlm.nih.gov/articles/PMC3327879/#:~:text=Schizophrenia%20is%20a%20severe%20mental,variation%20captured%20by%20common%20SNPs
            ),
        ),
        SumstatsSource(
            HAN_ASTHMA_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="asthma",
            pipe=CompositePipe(
                [
                    SetColToConstantPipe(
                        GWASLAB_SAMPLE_SIZE_COLUMN,
                        constant=393859,  # from summary statistics file,
                    ),  # True total sample size. From Gwas catalog
                    ComputeBetaPipe(),
                    ComputeSEPipe(),
                ]
            ),
            sample_info=BinaryPhenotypeSampleInfo(
                sample_prevalence=0.1638,  # 64,538 /(64,538 + 329,321)
                estimated_population_prevalence=0.117,  # See Johansson et al., 2019:  https://academic.oup.com/hmg/article/28/23/4022/5540983#:~:text=The%20disease%20prevalence%20in%20the%20Caucasian%20participants%20was%2011.7%25
            ),
        ),
        # https://www.nature.com/articles/s41467-020-15649-3#Sec28:~:text=The%20GWAS%20analysis%20in%20the%20UK%20Biobank%20included%2064%2C538%20cases%20and%20329%2C321%20controls
        SumstatsSource(
            KEATON_DBP_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="dbp",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                729882,  # from summary statistics file
            ),
            sample_info=QuantPhenotype(),
        ),
    ],
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    build="19",
)
