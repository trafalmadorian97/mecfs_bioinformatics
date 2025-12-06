from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.ldl_standard_sldsc import (
    LDL_STANDARD_SLDSC_TASK_GROUP,
)
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
        [
            MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.inner.bar_plot_task,
            MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.inner.labeled_filtered_gene_analysis_task,
            # MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC,
            # MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC_FDR_FILTERED,
            # LDL_BASE_CELL_PLOT,
            # PARTITIONED_MODEL_CAHOY_LD_SCORES_EXTRACTED,
            # PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_EXTRACTED,
            # PARTITIONED_MODEL_IMMGEN_LD_SCORES_EXTRACTED,
            # PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED,
        ]
        + LDL_STANDARD_SLDSC_TASK_GROUP.get_terminal_tasks(),
        incremental_save=True,
    )


if __name__ == "__main__":
    run_ldl_analysis()
