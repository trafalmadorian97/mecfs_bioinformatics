from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me_minus_pain.analysis.standard_analysis_decodeme_minus_pain_ols import (
    DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS,
)


def run_decode_me_minus_pain_ols():
    """
    Analysis of the residual after Johnston's GWAS of multisite pan is subtracted from
    DecodeME
    """
    DEFAULT_RUNNER.run(
        # DECODE_ME_MINUS_PAIN_OLS_GENETIC_CORR_GENERATOR.terminal_tasks()
        # + DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS.get_terminal_tasks()
        # [
            # DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS.magma_gene_manhattan_plot_unwrap
            DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS.h_magma_tasks.terminal_tasks()
        ,
        must_rebuild_transitive=[
            # DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.magma_independent_cluster_plot_unwrap.stepwise_cluster_list_task
        ],
    )


if __name__ == "__main__":
    run_decode_me_minus_pain_ols()
