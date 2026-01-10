from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc2022_sch_magma_human_brain_atlas_results_multiple_testing import (
    MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_MULTIPLE_TESTING,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc_2022_sch_magma_human_brain_atlas_forward_stepwise import (
    MAGMA_PGC2022_HBA_FORWARD_STEPWISE_CLUSTER_SELECT,
)
from mecfs_bio.build_system.task.magma.magma_plot_brain_atlas_result_with_stepwise_labels import (
    MAGMAPlotBrainAtlasResultWithStepwiseLabels,
)

PGC2022_HBA_MAGMA_STEPWISE_PLOT = MAGMAPlotBrainAtlasResultWithStepwiseLabels.create(
    result_table_task=MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_MULTIPLE_TESTING,
    asset_id="pgc2022_sch_human_brain_atlas_magma_stepwise_plot",
    stepwise_cluster_list_task=MAGMA_PGC2022_HBA_FORWARD_STEPWISE_CLUSTER_SELECT,
)
