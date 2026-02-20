"""
This file contains a basic driver script to analyze GWAS of the Alzheimer's phenotype.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.analysis.bellenguez_standard_analysis import (
    BELLENGUEZ_STANDARD_ANALYSIS,
)


def run_initial_alzheimers_analysis():
    """
    Function to analyze GWAS of Alzheimer's phenotype.
    Includes:
    - S-LDSC analysis
    - MAGMA analysis (Using both GTEx and HBA reference data)
    """
    DEFAULT_RUNNER.run(
        BELLENGUEZ_STANDARD_ANALYSIS.get_terminal_tasks()
        + [BELLENGUEZ_STANDARD_ANALYSIS.hba_magma_tasks.magma_independent_cluster_plot],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_initial_alzheimers_analysis()
