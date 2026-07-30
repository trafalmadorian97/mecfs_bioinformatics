from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.migraine.uk_biobank_2025.analysis.uk_biobank_2025_migraine_standard_analysis import (
    UK_BIOBANK_2025_EUR_MIGRAINE_STANDARD_ANALYSIS,
)


def run_initial_ukb_2025_migraine_analysis():
    """
    Function to run initial analysis on Million Veterans CFS phenotyp.

    Includes:

    - Generation of manhattan plot.
    - Extraction of lead variants.
    - Application of MAGMA to GTEx data to identify possible key tissues.
    - Application of MAGMA to human brain atlas to identify key brain tissue
    - Use of S-LDSC to investigate key cell and tissue types


    """
    DEFAULT_RUNNER.run(
        UK_BIOBANK_2025_EUR_MIGRAINE_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_initial_ukb_2025_migraine_analysis()
