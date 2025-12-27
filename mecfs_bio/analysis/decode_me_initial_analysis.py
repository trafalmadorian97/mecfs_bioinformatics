"""
Script to run initial analysis on DecodeME data.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_lead_variants_gene_labeled import (
    DECODE_ME_GWAS_1_GWASLAB_FILTERED_GENE_LIST_WITH_GENE_METADATA,
)


def run_initial_decode_me_analysis():
    """
    Function to run initial analysis on DecodeME data.

    Includes:

    - Generation of manhattan plot.
    - Extraction of lead variants.
    - Application of MAGMA to GTEx data to identify possible key tissues.
    - Use of S-LDSC to investigate key cell and tissue tpyes
    """
    DEFAULT_RUNNER.run(
        [
            # DECODE_ME_GWAS_1_MANHATTAN_AND_QQ_PLOT,
            # DECODE_ME_GWAS_1_MANHATTAN_PLOT,
            # DECODE_ME_GWAS_1_LEAD_VARIANTS,
            # DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST,
            # MAGMA_DECODE_ME_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
            DECODE_ME_GWAS_1_GWASLAB_FILTERED_GENE_LIST_WITH_GENE_METADATA
            # # DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA,
            # DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA_DROP_COLS,
        ],
        # + DECODE_ME_S_LDSC.get_terminal_tasks(),
        incremental_save=True,
    )


if __name__ == "__main__":
    run_initial_decode_me_analysis()
