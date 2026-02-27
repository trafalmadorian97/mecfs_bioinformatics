
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.analysis.hba_magma_analysis import (
    LEE_ET_AL_2018_HBA_MAGMA_TASKS_EDU,
)
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.analysis.lee_et_al_edu_cis_pqtl_mr import \
    LEE_ET_AL_EDU_CIS_PQTL_MR
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.analysis.lee_et_al_standard_s_ldsc import \
    LEE_ET_AL_EDU_STANDARD_SLDSC_TASK_GROUP
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import (
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS,
)


def run_edu_analysis():
    DEFAULT_RUNNER.run(
    LEE_ET_AL_EDU_CIS_PQTL_MR.terminal_tasks(),
        incremental_save=True,
    )


if __name__ == "__main__":
    run_edu_analysis()
