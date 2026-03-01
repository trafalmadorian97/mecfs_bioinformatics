from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_standard_analysis import (
    BENTHAM_LUPUS_STANDARD_ANALYSIS,
)


def lupus_analysis():
    """
    Script to analyze the lupus GWAS of Bentham et al
    """
    DEFAULT_RUNNER.run(
        [
            # BENTHAM_LUPUS_STANDARD_ANALYSIS.hba_magma_tasks.extracted_plot_task,
            # BENTHAM_LUPUS_STANDARD_ANALYSIS.hba_magma_tasks.magma_independent_cluster_plot
        ],
        BENTHAM_LUPUS_STANDARD_ANALYSIS.get_terminal_tasks(),
        # [BENTHAM_2015_RAW_BUILD_37] ,
        incremental_save=True,
        must_rebuild_transitive=[
            # BENTHAM_LUPUS_STANDARD_ANALYSIS.hba_magma_tasks.magma_independent_cluster_plot
            # BENTHAM_LUPUS_STANDARD_ANALYSIS.hba_magma_tasks.magma_hba_result_plot_task,
            # BENTHAM_LUPUS_STANDARD_ANALYSIS.hba_magma_tasks.magma_independent_cluster_plot
        ],
    )


if __name__ == "__main__":
    lupus_analysis()
