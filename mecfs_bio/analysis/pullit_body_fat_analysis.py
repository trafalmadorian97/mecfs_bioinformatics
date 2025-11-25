from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.body_fat_distribution.pullit_et_al_2018.processed_gwas_data.pullit_et_al_magma_task_generator import \
    PULLIT_ET_AL_2018_COMBINED_MAGMA_TASKS


"""
Plan: add a preprocessing pipe to get RSIDS that are embedded in SNPID
"""
def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [PULLIT_ET_AL_2018_COMBINED_MAGMA_TASKS.inner.bar_plot_task],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_initial_analysis()
