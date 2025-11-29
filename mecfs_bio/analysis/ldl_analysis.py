from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ldl.million_veterans.processed_gwas_data.million_veterans_ldl_eur_magma_task_generator import (
    MILLION_VETERANS_EUR_LDL_MAGMA_TASKS,
)


def run_ldl_analysis():
    """
    This script runs a basic analysis of the LDL GWAS from the Million Veterans Program

    This includes:

    - Using MAGMA to combine the LDL GWAS summary statistics with tissue-specific expression data from GTEx to identify key tissues involved in control of LDL levels.
    """
    DEFAULT_RUNNER.run(
        [MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.inner.bar_plot_task],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_ldl_analysis()
