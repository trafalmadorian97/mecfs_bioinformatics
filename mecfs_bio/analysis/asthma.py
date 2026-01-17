"""
Script to analyze asthma phenotype.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.mr.asthma_ppp_cis_mr_markdown import (
    ASTHMA_HAN_PPP_CIS_MR_MARKDOWN,
)


def run_initial_asthma_analysis():
    """
    Function to analyze asthma phenotype
    Includes:
    - MAGMA and S-LDSC analysis
    """
    DEFAULT_RUNNER.run(
        [
            # ASTHMA_PP_CIS_MR_UNIPROT_LABELED,
            ASTHMA_HAN_PPP_CIS_MR_MARKDOWN
        ],
        # [HAN_2022_ASTHMA_TSMR],
        # [HAN_ET_AL_ASTHMA_RAW] + HAN_ASTHMA_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            ASTHMA_HAN_PPP_CIS_MR_MARKDOWN
            # ASTHMA_PP_CIS_MR_UNIPROT_LABELED
            # HAN_ASTHMA_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task
            # HAN_2022_ASTHMA_TSMR
        ],
    )


if __name__ == "__main__":
    run_initial_asthma_analysis()
