"""
Analogous to MAGMA_DECODE_ME_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT, but controls for average across brain tissue, not all tissues.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_gwas_1_build_37_magma_ensembl_specific_tissue_gene_sets_brain_average import (
    MAGMA_DECODE_ME_SPECIFIC_TISSUE_GENE_SET_ANALYSIS_BRAIN_AVERAGE,
)
from mecfs_bio.build_system.task.magma.magma_plot_gene_set_result import (
    MAGMAPlotGeneSetResult,
)

MAGMA_DECODE_ME_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT_BRAIN_AVERAGE = MAGMAPlotGeneSetResult.create(
    gene_set_analysis_task=MAGMA_DECODE_ME_SPECIFIC_TISSUE_GENE_SET_ANALYSIS_BRAIN_AVERAGE,
    asset_id="decode_me_gwas_1_build_37_magma_ensemble_specific_tissue_gene_covar_analysis_bar_plot_brain_average",
)


