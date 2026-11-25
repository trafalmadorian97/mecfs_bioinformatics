"""
Script to analyze the device-measured sleep GWAS of Portas et al. (2026).
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.sleep.portas_et_al_2026.analysis.portas_sleep_standard_analysis import (
    PORTAS_NREM_SLEEP_DURATION_STANDARD_ANALYSIS,
    PORTAS_REM_SLEEP_DURATION_STANDARD_ANALYSIS,
    PORTAS_SLEEP_DURATION_STANDARD_ANALYSIS,
    PORTAS_SLEEP_EFFICIENCY_STANDARD_ANALYSIS,
)


def run_portas_sleep_analysis() -> None:
    """
    Apply the standard analysis pipeline to the sleep duration, sleep efficiency, REM sleep duration,
    and non-REM sleep duration GWAS of Portas et al. (2026).
    Analysis includes:
    - Extraction of lead variants.
    - Application of MAGMA to GTEx data to identify possible key tissues.
    - Application of MAGMA to human brain atlas to identify key brain tissue
    - Use of S-LDSC to investigate key cell and tissue types
    - Estimation of SNP heritability via LDSC
    """
    DEFAULT_RUNNER.run(
        PORTAS_SLEEP_DURATION_STANDARD_ANALYSIS.get_terminal_tasks()
        + PORTAS_SLEEP_EFFICIENCY_STANDARD_ANALYSIS.get_terminal_tasks()
        + PORTAS_REM_SLEEP_DURATION_STANDARD_ANALYSIS.get_terminal_tasks()
        + PORTAS_NREM_SLEEP_DURATION_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_portas_sleep_analysis()
