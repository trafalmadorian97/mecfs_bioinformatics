from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.analysis.magma.liu_et_al_2023_eur_37_specific_tissue_bar_plot import (
    LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
)


def run_ibd_analysis():
    """
    This script runs a basic analysis of the IBD GWAS from Liu et al.
    It includes:


    - Using MAGMA to combine the IBD GWAS summary statistics with tissue-specific expression data from GTEx to identify possible tissues involved in the IBD disease process.
    """
    DEFAULT_RUNNER.run(
        [LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT], incremental_save=True
    )


if __name__ == "__main__":
    run_ibd_analysis()
