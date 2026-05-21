from mecfs_bio.asset_generator.genetic_correlation_asset_generator import (
    genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.million_veterans_ldl_eur_magma_task_generator import (
    MILLION_VETERANS_EUR_LDL_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.ldl.million_veterans.auxiliary.mv_ldl_pheotype_info import (
    MV_LDL_PHENOTYPE_INFO,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_standard_analysis import (
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.auxiliary.willer_ldl_phenotype_info import (
    WILLER_LDL_PHENOTYPE_INFO,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    SumstatsSource,
)

CT_LDSC_LDL = genetic_corr_by_ct_ldsc_asset_generator(
    base_name="ldl_ct_ldsc",
    sources=[
        SumstatsSource(
            MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.sumstats_task,
            alias="Million_Veterans_LDL",
            sample_info=MV_LDL_PHENOTYPE_INFO,
        ),
        SumstatsSource(
            WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.magma_tasks.sumstats_task,
            alias="Willer_et_Al_LDL",
            sample_info=WILLER_LDL_PHENOTYPE_INFO,
        ),
    ],
)
