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
        BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks["immgen"].plot_task,
            BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks["multi_tissue_chromatin"].plot_task,
            BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks["multi_tissue_gene_expression"].plot_task,
        ],
        incremental_save=True,
        must_rebuild_transitive=[
            BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks["immgen"].plot_task,
            BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks["multi_tissue_chromatin"].plot_task,
            BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks["multi_tissue_gene_expression"].plot_task,
        ],
    )


if __name__ == "__main__":
    lupus_analysis()
