"""
This file contains a basic driver script to analyze the GWAS of the multisite pain phenotype.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)


def run_initial_multisite_pain_analysis():
    """
    Function to analyze GWAS of multisite pain phenotype.
    Includes:
    - S-LDSC analysis
    - MAGMA analysis (Using both GTEx and HBA reference data)
    """
    DEFAULT_RUNNER.run(
        JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
    )


if __name__ == "__main__":
    run_initial_multisite_pain_analysis()
