"""
Script to analyze inflammatory bowl disease GWAS from Liu et al.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.analysis.magma.liu_et_al_2023_eur_37_hba_magma_analysis import (
    IBD_HBA_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.analysis.magma.liu_et_al_2023_eur_37_specific_tissue_bar_plot import (
    LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
)


def run_ibd_analysis():
    """
    This script runs a basic analysis of the IBD GWAS from Liu et al.
    It includes:


    - Using MAGMA to combine the IBD GWAS summary statistics with tissue-specific expression data from GTEx to identify possible tissues involved in the IBD disease process.
    - S-LDSC analysis using standard reference data
    """
    DEFAULT_RUNNER.run(
        [
            LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
        ]
        + IBD_HBA_MAGMA_TASKS.terminal_tasks(),
        # + LIU_ET_AL_S_LSDC_FROM_SNP_150.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[IBD_HBA_MAGMA_TASKS.magma_independent_cluster_plot],
    )


if __name__ == "__main__":
    run_ibd_analysis()
