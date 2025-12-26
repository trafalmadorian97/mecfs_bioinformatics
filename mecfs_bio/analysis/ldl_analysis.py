"""
Script to run basic analysis of LDL GWAS from the Million Veterans Program.
"""
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.ldl_standard_sldsc import (
    LDL_STANDARD_SLDSC_TASK_GROUP,
)
from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.million_veterans_ldl_eur_magma_task_generator import (
    MILLION_VETERANS_EUR_LDL_MAGMA_TASKS,
)


def run_ldl_analysis():
    """
    This script runs a basic analysis of the LDL GWAS from the Million Veterans Program

    This includes:

    - Using MAGMA to combine the LDL GWAS summary statistics with tissue-specific expression data from GTEx to identify key tissues involved in control of LDL levels.
    - Using S-LDSC to identify key tissues and cell types via a heritability model
    """
    DEFAULT_RUNNER.run(
        [
            MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.inner.bar_plot_task,
            MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.inner.labeled_filtered_gene_analysis_task,
        ]
        + LDL_STANDARD_SLDSC_TASK_GROUP.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            LDL_STANDARD_SLDSC_TASK_GROUP.partitioned_tasks[
                "multi_tissue_gene_expression"
            ].add_categories_task
        ],
    )


if __name__ == "__main__":
    run_ldl_analysis()
