"""
Script to analyze imaging-derived right-heart phenotypes
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.imaging_derived_heart_phenotypes.pirruccello_et_al_2022.analysis.rvef_standard_analysis import (
    RVEF_STANDARD_ANALYSIS_ASSIGN_RSID,
)


def run_initial_right_heart_analysis():
    """
    Function to analyze imaging-derived right-heart phenotype.
    Includes:
    - MAGMA and S-LDSC analysis of LVEF
    """
    DEFAULT_RUNNER.run(
        RVEF_STANDARD_ANALYSIS_ASSIGN_RSID.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_initial_right_heart_analysis()
