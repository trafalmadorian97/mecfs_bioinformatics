from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.analysis.verweiji_standard_analysis import (
    VERWEIJI_ET_AL_HRR_STANDARD_ANALYSIS,
)


def run_heart_rate_recovery_analysis():
    """
    Analyze GWAS of heart rate recovery
    uses GWAS summary statistics from Verweiji et al. (2018)

    """
    DEFAULT_RUNNER.run(
        VERWEIJI_ET_AL_HRR_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_heart_rate_recovery_analysis()
