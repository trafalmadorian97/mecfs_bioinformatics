"""
Script to analyze the body fat distribution GWAS from Pullit et al.

Full citation:
Pulit, Sara L., et al. "Meta-analysis of genome-wide association studies for body fat distribution in 694 649 individuals of European ancestry." Human molecular genetics 28.1 (2019): 166-174.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.body_fat_distribution.pullit_et_al_2018.processed_gwas_data.pullit_et_al_magma_task_generator import (
    PULLIT_ET_AL_2018_COMBINED_MAGMA_TASKS,
)


def run_body_fat_distribution_analysis():
    """
    Analyze body fat distribution GWAS from pullit et al.
    Includes:
    - MAGMA analysis of key tissues using GTEx bulk RNAseq reference data

    """
    DEFAULT_RUNNER.run(
        [PULLIT_ET_AL_2018_COMBINED_MAGMA_TASKS.inner.bar_plot_task],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_body_fat_distribution_analysis()
