"""
Script to run initial analysis on Lee et al. Educational Attainment GWAS
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.analysis.hba_magma_analysis import (
    LEE_ET_AL_2018_HBA_MAGMA_TASKS_EDU,
)
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import (
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS,
)


def run_edu_analysis():
    """
    Function to run initial analysis on the Lee et al. Educational Attainment GWAS
    Includes:

    - Calculating heritability by LDSC
    - Plotting key tissues as determined by MAGMA Gene Level Analysis of key cells and tissues of SLDSC
    - Use S-LSDC to look for key tissues and cells using standard reference datasets

    """
    DEFAULT_RUNNER.run(
        [
            LEE_ET_AL_2018_COMBINED_MAGMA_TASKS.inner.bar_plot_task,
            # LEE_ET_AL_2018_H2_BY_LDSC,
        ]
        + LEE_ET_AL_2018_HBA_MAGMA_TASKS_EDU.terminal_tasks(),
        # + LEE_ET_AL_EDU_STANDARD_SLDSC_TASK_GROUP.get_terminal_tasks(),
        incremental_save=True,
    )


if __name__ == "__main__":
    run_edu_analysis()
