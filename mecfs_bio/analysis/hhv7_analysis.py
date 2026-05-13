"""
Script to perform standard analysis on Kamitaki et al.'s GWAS of HHV7 DNA levels.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.human_herpesvirus_7_dna.kamitaki_et_al_2025.analysis.kamitaki_et_al_2025_standard_analysis import (
    KAMITAKI_ET_AL_STANDARD_ANALYSIS,
)


def run_hhv7():
    """
    Function to perform standard analysis on Kamitaki et al.'s GWAS of HHV7 DNA levels.

    Includes S-LDSC, GTEx MAGMA, and HBA MAGMA
    """
    DEFAULT_RUNNER.run(
        KAMITAKI_ET_AL_STANDARD_ANALYSIS.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_hhv7()
