from mecfs_bio.asset_generator.genetic_correlation_asset_generator import (
    genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13.analysis.finngen_r13_migraine_standard_analysis import (
    FINNGEN_R13_MIGRAINE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13.auxiliary.prevalence_info import (
    FINNGEN_R13_MIGRAINE_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_no_aura.analysis.finngen_r13_migraine_no_aura_standard_analysis import (
    FINNGEN_R13_MIGRAINE_NO_AURA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_no_aura.auxiliary.prevalence_info import (
    FINNGEN_R13_MIGRAINE_NO_AURA_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_triptan.analysis.finngen_r13_migraine_triptan_standard_analysis import (
    FINNGEN_R13_MIGRAINE_TRIPTAN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_triptan.auxiliary.prevalence_info import (
    FINNGEN_R13_MIGRAINE_TRIPTAN_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_with_aura.analysis.finngen_r13_migraine_with_aura_standard_analysis import (
    FINNGEN_R13_MIGRAINE_WITH_AURA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_with_aura.auxiliary.prevalence_info import (
    FINNGEN_R13_MIGRAINE_WITH_AURA_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.migraine.million_veterans.analysis.million_veterans_migraine_standard_analysis import (
    MILLION_VETERANS_EUR_MIGRAINE_STANDARD_ANALYSIS,
    MILLION_VETERANS_MIGRAINE_SAMPLE_INFO,
)
from mecfs_bio.assets.gwas.migraine.uk_biobank_2025.analysis.uk_biobank_2025_migraine_standard_analysis import (
    UK_BIOBANK_2025_EUR_MIGRAINE_STANDARD_ANALYSIS,
    UK_BIOBANK_MIGRAINE_SAMPLE_INFO,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    SumstatsSource,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe

MIGRAINE_CROSS_STUDY_CT_LDSC_ASSET_GENERATOR = genetic_corr_by_ct_ldsc_asset_generator(
    base_name="migraine_ct_ldsc",
    sources=[
        SumstatsSource(
            UK_BIOBANK_2025_EUR_MIGRAINE_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="UK_Biobank",
            sample_info=UK_BIOBANK_MIGRAINE_SAMPLE_INFO,
        ),
        SumstatsSource(
            MILLION_VETERANS_EUR_MIGRAINE_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="Million_Veterans",
            pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
            sample_info=MILLION_VETERANS_MIGRAINE_SAMPLE_INFO,
        ),
        SumstatsSource(
            FINNGEN_R13_MIGRAINE_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="FinnGen_Migraine",
            sample_info=FINNGEN_R13_MIGRAINE_PREVALENCE_INFO,
        ),
        SumstatsSource(
            FINNGEN_R13_MIGRAINE_TRIPTAN_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="FinnGen_Migraine_Triptan",
            sample_info=FINNGEN_R13_MIGRAINE_TRIPTAN_PREVALENCE_INFO,
        ),
        SumstatsSource(
            FINNGEN_R13_MIGRAINE_NO_AURA_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="FinnGen_Migraine_No_Aura",
            sample_info=FINNGEN_R13_MIGRAINE_NO_AURA_PREVALENCE_INFO,
        ),
        SumstatsSource(
            FINNGEN_R13_MIGRAINE_WITH_AURA_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="FinnGen_Migraine_With_Aura",
            sample_info=FINNGEN_R13_MIGRAINE_WITH_AURA_PREVALENCE_INFO,
        ),
    ],
)
