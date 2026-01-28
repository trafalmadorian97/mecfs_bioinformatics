"""
Script to run initial analysis on DecodeME data.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_combined_gene_list_markdown import (
    DECODE_ME_MASTER_GENE_LIST_AS_MARKDOWN,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_lead_variants import (
    DECODE_ME_GWAS_1_LEAD_VARIANTS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_manhattan import (
    DECODE_ME_GWAS_1_MANHATTAN_PLOT,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_manhattan_and_qq import (
    DECODE_ME_GWAS_1_MANHATTAN_AND_QQ_PLOT,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import (
    DECODE_ME_S_LDSC,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_hba_magma_analysis import (
    DECODE_ME_HBA_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.magma_specific_tissue_bar_plot import (
    MAGMA_DECODE_ME_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_cis_pqtl_mr import (
    DECODE_ME_BASIC_CIS_PQTL_MR,
)


def run_initial_decode_me_analysis():
    """
    Function to run initial analysis on DecodeME data.

    Includes:

    - Generation of manhattan plot.
    - Extraction of lead variants.
    - Application of MAGMA to GTEx data to identify possible key tissues.
    - Use of S-LDSC to investigate key cell and tissue types
    """
    DEFAULT_RUNNER.run(
        # [
        DECODE_ME_GWAS_1_MANHATTAN_AND_QQ_PLOT,
        DECODE_ME_GWAS_1_MANHATTAN_PLOT,
        DECODE_ME_GWAS_1_LEAD_VARIANTS,
        #     DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST,
        MAGMA_DECODE_ME_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
        #     # DECODE_ME_GWAS_1_COMBINED_GENE_LISTS,
        DECODE_ME_MASTER_GENE_LIST_AS_MARKDOWN,
        #     # # DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA,
        #     # DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA_DROP_COLS,
        #
        # ]
        +DECODE_ME_S_LDSC.get_terminal_tasks()
        + DECODE_ME_HBA_MAGMA_TASKS.terminal_tasks()
        + DECODE_ME_BASIC_CIS_PQTL_MR.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            # DECODE_ME_MASTER_GENE_LIST_AS_MARKDOWN,
            # DECODE_ME_MASTER_GENE_LIST_WITH_GGET
        ],
    )


if __name__ == "__main__":
    run_initial_decode_me_analysis()
