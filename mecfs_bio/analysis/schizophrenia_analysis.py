"""
Analyze Schizophrenia dataset from PGC 2022
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import (
    SCH_PGC_2022_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.raw.raw_sch_pgc2022 import (
    PGC_2022_SCH_RAW,
)


def run_initial_schizophrenia_analysis():
    """
    Function to analyze GWAS of Schizophrenia phenotype.
    Includes:
    - MAGMA and S-LDSC analysis
    """
    DEFAULT_RUNNER.run(
        [
            PGC_2022_SCH_RAW,
        ]
        + SCH_PGC_2022_STANDARD_ANALYSIS.get_terminal_tasks()
    )


if __name__ == "__main__":
    run_initial_schizophrenia_analysis()
