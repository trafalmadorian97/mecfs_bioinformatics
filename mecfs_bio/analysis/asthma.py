"""
Script to analyze asthma phenotype.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER


def run_initial_asthma_analysis():
    """
    Function to analyze asthma phenotype
    Includes:
    - MAGMA and S-LDSC analysis
    """
    DEFAULT_RUNNER.run(
        [],
        # [HAN_2022_ASTHMA_TSMR],
        # [HAN_ET_AL_ASTHMA_RAW] + HAN_ASTHMA_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            # HAN_ASTHMA_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task
            # HAN_2022_ASTHMA_TSMR
        ],
    )


if __name__ == "__main__":
    run_initial_asthma_analysis()
