"""
Script to analyze heart rate recovery GWAS of Verweiji et al.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.analysis.verweiji_standard_analysis import (
    VERWEIJI_ET_AL_HRR_STANDARD_ANALYSIS,
)


def run_heart_rate_recovery_analysis():
    """
    Analyze GWAS of heart rate recovery ( Verweiji et al. (2018))
    Analysis includes:
    - Application of MAGMA to GTEx data to identify possible key tissues.
    - Application of MAGMA to human brain atlas to identify key brain tissue
    - Use of S-LDSC to investigate key cell and tissue types

    """
    DEFAULT_RUNNER.run(
        VERWEIJI_ET_AL_HRR_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_heart_rate_recovery_analysis()
