"""
Analyze Schizophrenia dataset from PGC 2022
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc_2022_sch_human_brain_atlas_stepwise_plot import (
    PGC2022_HBA_MAGMA_STEPWISE_PLOT,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc_2022_sch_magma_human_brain_atlas_forward_stepwise import (
    MAGMA_PGC2022_HBA_FORWARD_STEPWISE_CLUSTER_SELECT,
)


def run_initial_schizophrenia_analysis():
    """
    Function to analyze GWAS of Schizophrenia phenotype.
    Includes:
    - S-LDSC analysis
    - MAGMA analysis (Using both GTEx and HBA reference data)
    """
    DEFAULT_RUNNER.run(
        [
            # PGC_2022_SCH_RAW,
            # PGC2022_SCH_MAGMA_ENTREZ_ANNOTATIONS,
            # PGC2022_SCH_MAGMA_ENTREZ_GENE_ANALYSIS
            # MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR
            # PGC_2022_SCH_MAGMA_HUMAN_BRAIN_ATLAS_PLOT,
            # MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_CONDITIONAL_ANALYSIS
            # MAGMA_PGC2022_HBA_FORWARD_STEPWISE_CLUSTER_SELECT,
            PGC2022_HBA_MAGMA_STEPWISE_PLOT
        ],
        # + SCH_PGC_2022_STANDARD_ANALYSIS.get_terminal_tasks(),
        # must_rebuild_transitive=[PGC2022_SCH_MAGMA_ENTREZ_ANNOTATIONS]
        # must_rebuild_transitive=[PGC2022_HBA_MAGMA_SPEC_MATRIX_FILTERED],
        must_rebuild_transitive=[
            PGC2022_HBA_MAGMA_STEPWISE_PLOT,
            MAGMA_PGC2022_HBA_FORWARD_STEPWISE_CLUSTER_SELECT,
        ],
    )


if __name__ == "__main__":
    run_initial_schizophrenia_analysis()
