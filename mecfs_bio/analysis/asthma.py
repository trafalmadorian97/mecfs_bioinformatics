from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.asthma_standard_analysis import (
    HAN_ASTHMA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.asthma.han_et_al_2022.raw.raw_asthma_data import (
    HAN_ET_AL_ASTHMA_RAW,
)


def run_initial_asthma_analysis():
    """
    Function to analyze asthma phenotype
    Includes:
    - MAGMA and S-LDSC analysis
    """
    DEFAULT_RUNNER.run(
        [HAN_ET_AL_ASTHMA_RAW] + HAN_ASTHMA_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            HAN_ASTHMA_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task
        ]
    )


if __name__ == "__main__":
    run_initial_asthma_analysis()
