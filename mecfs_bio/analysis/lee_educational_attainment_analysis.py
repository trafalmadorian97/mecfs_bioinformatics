from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.analysis.lee_et_al_h2_by_ldsc import (
    LEE_ET_AL_2018_H2_BY_LDSC,
)
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.processed_gwas_data.lee_et_al_magma_task_generator import (
    LEE_ET_AL_2018_COMBINED_MAGMA_TASKS,
)


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [
            LEE_ET_AL_2018_COMBINED_MAGMA_TASKS.inner.bar_plot_task,
            LEE_ET_AL_2018_H2_BY_LDSC,
        ]
    )


if __name__ == "__main__":
    run_initial_analysis()
