"""
Task generator to use cross-trait LD Score Regression to estimate genetic correlation between a several traits.
Each pair of traits is a Task.  This facilitates caching of results
"""

from mecfs_bio.asset_generator.genetic_correlation_asset_generator import (
    genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.asthma_standard_analysis import (
    HAN_ASTHMA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.asthma.han_et_al_2022.auxiliary.prevalence_info import (
    HAN_ET_AL_ASTHMA_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import (
    KEATON_DBP_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import (
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.analysis.verweiji_standard_analysis import (
    VERWEIJI_ET_AL_HRR_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.auxiliary.prevalence_nfo import (
    LIU_ET_AL_IBD_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_with_rsid_sumstats import (
    LIU_ET_AL_SUMSTATS_WITH_RSID_FROM_SNP150,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_with_annovar_37_rsids_sumstats import (
    DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.auxiliary.prevalance_info import (
    PGC_2022_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import (
    SCH_PGC_2022_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.auxiliary.prevelnce_info import (
    AEGISDOTTIR_SYNCOPE_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_sumstats_liftover_hapmap3_dedup import (
    AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS_HAPMAP3_DEDUP,
)
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_standard_analysis import (
    BENTHAM_LUPUS_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.auxiliary.prevalence_info import (
    BENTHAM_LUPUS_PREVALENCE_INFO,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
    SumstatsSource,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import (
    ComputeBetaIfNeededPipe,
    ComputeBetaPipe,
)
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

CT_LDSC_INITIAL_ASSET_GENERATOR = genetic_corr_by_ct_ldsc_asset_generator(
    "initial_rg",
    sources=[
        SumstatsSource(
            DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
            alias="DecodeME",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                constant=275488,
            ),  # true total sample size. From preprint
            sample_info=DECODE_ME_PREVALENCE_INFO,
        ),
        SumstatsSource(
            VERWEIJI_ET_AL_HRR_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="HR_recovery",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                58_818,  # from abstract
            ),
            sample_info=QuantPhenotype(),
        ),
        SumstatsSource(
            JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="Multisite_pain",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN, constant=387649
            ),  # True total sample size. From Gwas catalog
            sample_info=QuantPhenotype(),
        ),
        SumstatsSource(
            SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="Schizophrenia",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                constant=130644,  # from summary statistics file (130644=53386+77258)
            ),  # True total sample size. From Gwas catalog
            sample_info=PGC_2022_PREVALENCE_INFO,
        ),
        SumstatsSource(
            HAN_ASTHMA_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="Asthma",
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
            sample_info=HAN_ET_AL_ASTHMA_PREVALENCE_INFO,
        ),
        # https://www.nature.com/articles/s41467-020-15649-3#Sec28:~:text=The%20GWAS%20analysis%20in%20the%20UK%20Biobank%20included%2064%2C538%20cases%20and%20329%2C321%20controls
        SumstatsSource(
            KEATON_DBP_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="DBP",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                729882,  # from summary statistics file
            ),
            sample_info=QuantPhenotype(),
        ),
        SumstatsSource(
            LEE_ET_AL_2018_COMBINED_MAGMA_TASKS.sumstats_task,
            alias="Educational_attainment",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                257841,  # source: GWASLAB metadata
            ),
            sample_info=QuantPhenotype(),
        ),
        SumstatsSource(
            LIU_ET_AL_SUMSTATS_WITH_RSID_FROM_SNP150,
            alias="Inflammatory_bowel_disease",
            pipe=SetColToConstantPipe(
                GWASLAB_SAMPLE_SIZE_COLUMN,
                59957,  # 25,042 + 34,915
            ),
            sample_info=LIU_ET_AL_IBD_PREVALENCE_INFO,
        ),
        SumstatsSource(
            AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS_HAPMAP3_DEDUP,
            alias="Syncope",
            pipe=CompositePipe(
                [
                    SetColToConstantPipe(
                        GWASLAB_SAMPLE_SIZE_COLUMN,
                        946_861,
                    ),
                    ComputeBetaIfNeededPipe(),
                ]
            ),
            sample_info=AEGISDOTTIR_SYNCOPE_PREVALENCE_INFO,
        ),
        SumstatsSource(
            BENTHAM_LUPUS_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="Lupus",
            pipe=CompositePipe(
                [SetColToConstantPipe(GWASLAB_SAMPLE_SIZE_COLUMN, 14267)]
            ),
            sample_info=BENTHAM_LUPUS_PREVALENCE_INFO,
        ),
    ],
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    build="19",
)
