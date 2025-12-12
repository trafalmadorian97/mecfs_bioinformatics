from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.analysis.magma.liu_et_al_2023_eur_37_magma_ensembl_specific_tisssue_gene_sets import (
    LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_GENE_SET_ANALYSIS,
)
from mecfs_bio.build_system.task.magma.magma_plot_gene_set_result import (
    MAGMAPlotGeneSetResult,
)

LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT = MAGMAPlotGeneSetResult.create(
    gene_set_analysis_task=LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_GENE_SET_ANALYSIS,
    asset_id="liu_et_Al_ibd_eur_build_37_magma_ensemble_specific_tissue_gene_covar_analysis_bar_plot",
)
