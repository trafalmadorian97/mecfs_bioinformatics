from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ebv_dna.nyeo_et_al.analysis.ebv_dna_standard_analysis import EBV_DNA_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.ebv_dna.nyeo_et_al.raw.ebv_dna import NYEO_EBV_DNA_SUMSTATS


def run_initial_ebv_analysis():
    DEFAULT_RUNNER.run(
        [EBV_DNA_STANDARD_ANALYSIS.tasks.hba_magma_tasks.extracted_plot_task],
        must_rebuild_transitive=[
            EBV_DNA_STANDARD_ANALYSIS.tasks.hba_magma_tasks.magma_hba_result_plot_task,
        ],
        incremental_save=True
    )


if __name__ == "__main__":
    run_initial_ebv_analysis()
