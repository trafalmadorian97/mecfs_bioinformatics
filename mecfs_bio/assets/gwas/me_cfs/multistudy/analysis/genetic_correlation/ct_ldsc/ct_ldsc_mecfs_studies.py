from mecfs_bio.asset_generator.genetic_correlation_asset_generator import (
    genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_with_annovar_37_rsids_sumstats import (
    DECODEME_ME_SUMSTATS_37_WITH_ANNOVAR_RSID,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_female.analysis.decode_me_female_standard_analysis import (
    DECODE_ME_FEMALE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_female.auxiliary.decode_me_female_phenotype_info import (
    DECODE_ME_FEMALE_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_gwas_2.analysis.decode_me_gwas_2_standard_analysis import (
    DECODE_ME_GWAS_2_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_gwas_2.auxiliary.decode_me_gwas_2_phenotype_info import (
    DECODE_ME_GWAS_2_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_infectious_onset.analysis.decode_me_infectious_onset_standard_analysis import (
    DECODE_ME_INFECTIOUS_ONSET_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_infectious_onset.auxiliary.decode_me_infectious_onset_phenotype_info import (
    DECODE_ME_INFECTIOUS_ONSET_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_male.analysis.decode_me_male_standard_analysis import (
    DECODE_ME_MALE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_male.auxiliary.decode_me_male_phenotype_info import (
    DECODE_ME_MALE_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_non_infectious_onset.analysis.decode_me_non_infectious_onset_standard_analysis import (
    DECODE_ME_NON_INFECTIOUS_ONSET_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_non_infectious_onset.auxiliary.decode_me_non_infectious_onset_phenotype_info import (
    DECODE_ME_NON_INFECTIOUS_ONSET_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.million_veterans.analysis.million_veterans_cfs_standard_analysis import (
    MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP,
)
from mecfs_bio.assets.gwas.me_cfs.million_veterans.auxiliary.pheotype_info import (
    MILLION_VETERAN_CFS_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.neale_lab.analysis.neale_lab_cfs_standard_analysis import (
    NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP,
)
from mecfs_bio.assets.gwas.me_cfs.neale_lab.auxiliary.neale_lab_cfs_phenotype_info import (
    NEALE_LAB_CFS_PHENOTYPE_INFO,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    SumstatsSource,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

CFS_CT_LDSC_ASSET_GENERATOR = genetic_corr_by_ct_ldsc_asset_generator(
    base_name="cfs_ct_ldsc",
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
            DECODE_ME_FEMALE_STANDARD_ANALYSIS.tasks.magma_tasks.sumstats_task,
            alias="DecodeME_Female",
            sample_info=DECODE_ME_FEMALE_PHENOTYPE_INFO,
        ),
        SumstatsSource(
            DECODE_ME_MALE_STANDARD_ANALYSIS.tasks.magma_tasks.sumstats_task,
            alias="DecodeME_Male",
            sample_info=DECODE_ME_MALE_PHENOTYPE_INFO,
        ),
        SumstatsSource(
            DECODE_ME_GWAS_2_STANDARD_ANALYSIS.tasks.magma_tasks.sumstats_task,
            alias="DecodeME_GWAS_2",
            sample_info=DECODE_ME_GWAS_2_PHENOTYPE_INFO,
        ),
        SumstatsSource(
            DECODE_ME_INFECTIOUS_ONSET_STANDARD_ANALYSIS.tasks.magma_tasks.sumstats_task,
            alias="DecodeME_Infectious_Onset",
            sample_info=DECODE_ME_INFECTIOUS_ONSET_PHENOTYPE_INFO,
        ),
        SumstatsSource(
            DECODE_ME_NON_INFECTIOUS_ONSET_STANDARD_ANALYSIS.tasks.magma_tasks.sumstats_task,
            alias="DecodeME_Non_Infectious_Onset",
            sample_info=DECODE_ME_NON_INFECTIOUS_ONSET_PHENOTYPE_INFO,
        ),
        SumstatsSource(
            MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP.magma_tasks.sumstats_task,
            alias="Million_Veterans_CFS",
            sample_info=MILLION_VETERAN_CFS_PHENOTYPE_INFO,
            pipe=CompositePipe([ComputeBetaIfNeededPipe(), ComputeSEPipe()]),
        ),
        SumstatsSource(
            NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP.magma_tasks.sumstats_task,
            alias="Neale_Lab_CFS",
            sample_info=NEALE_LAB_CFS_PHENOTYPE_INFO,
        ),
    ],
)
