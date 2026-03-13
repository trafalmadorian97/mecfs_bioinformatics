from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.magma.johnson_drg_magma_analysis_cepo import (
    MAGMA_JOHNSON_DRG_ANALYSIS_CEPO,
)
from mecfs_bio.build_system.task.magma.magma_plot_gene_set_result import (
    MAGMAPlotGeneSetResult,
)

JOHNSON_DRG_MAGMA_CEPO_BAR_PLOT = MAGMAPlotGeneSetResult.create(
    gene_set_analysis_task=MAGMA_JOHNSON_DRG_ANALYSIS_CEPO,
    asset_id="johnson_drg_magma_cepo_bar_plot",
    number_of_bars=15,
    label_col="VARIABLE",
)
