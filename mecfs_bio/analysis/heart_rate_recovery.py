from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.analysis.verweiji_lead_variant_tasks import (
    VERWEIJI_LEAD_VARIANTS,
)
from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.processed.verweiji_magma_task_generator import (
    VERWEIJI_ET_AL_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.raw.verweij_raw_summary import (
    VERWEIJI_ET_AL_RAW_HARMONIZED_BUILD_37,
)


def run_heart_rate_recovery_analysis():
    """
    Analyze GWAS of heart rate recovery
    uses GWAS summary statistics from Verweiji et al. (2018)

    """
    DEFAULT_RUNNER.run(
        [VERWEIJI_ET_AL_RAW_HARMONIZED_BUILD_37, VERWEIJI_LEAD_VARIANTS]
        + VERWEIJI_ET_AL_COMBINED_MAGMA_TASKS.inner.terminal_tasks(),
        # + VERWEIJI_SLDSC_TASK_GROUP.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            VERWEIJI_ET_AL_COMBINED_MAGMA_TASKS.inner.gget_labeled_filtered_gene_analysis_task
        ],
    )


if __name__ == "__main__":
    run_heart_rate_recovery_analysis()
