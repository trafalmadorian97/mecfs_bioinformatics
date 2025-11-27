from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.psoriasis.million_veterans.processed_gwas_data.million_veterans_psoriasis_eur_magma_task_generator import (
    MILLION_VETERANS_EUR_PSORIASIS_COMBINED_MAGMA_TASKS,
)


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [MILLION_VETERANS_EUR_PSORIASIS_COMBINED_MAGMA_TASKS.inner.bar_plot_task],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_initial_analysis()
