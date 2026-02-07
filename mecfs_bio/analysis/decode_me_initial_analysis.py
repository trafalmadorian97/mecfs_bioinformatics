"""
Script to run initial analysis on DecodeME data.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_manhattan import (
    DECODE_ME_GWAS_1_MANHATTAN_PLOT,
)


def run_initial_decode_me_analysis():
    """
    Function to run initial analysis on DecodeME data.

    Includes:

    - Generation of manhattan plot.
    - Extraction of lead variants.
    - Application of MAGMA to GTEx data to identify possible key tissues.
    - Application of MAGMA to human brain atlas to identify key brain tissue
    - Use of S-LDSC to investigate key cell and tissue types
    - USe of Mendelian Randomization to look for causal blood proteins


    Comment out lines to remove tasks
    """
    DEFAULT_RUNNER.run(
        [
            # DECODE_ME_GWAS_1_MANHATTAN_AND_QQ_PLOT,
            DECODE_ME_GWAS_1_MANHATTAN_PLOT,
            # DECODE_ME_GWAS_1_LEAD_VARIANTS,
            # MAGMA_DECODE_ME_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
            # DECODE_ME_MASTER_GENE_LIST_AS_MARKDOWN,
        ],
        # + DECODE_ME_S_LDSC.get_terminal_tasks()
        # + DECODE_ME_HBA_MAGMA_TASKS.terminal_tasks(),
        # + DECODE_ME_BASIC_CIS_PQTL_MR.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_initial_decode_me_analysis()
