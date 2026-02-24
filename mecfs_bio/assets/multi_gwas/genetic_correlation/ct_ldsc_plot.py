from mecfs_bio.assets.multi_gwas.genetic_correlation.ct_ldsc_initial import (
    CT_LDSC_INITIAL,
)
from mecfs_bio.build_system.task.genetic_correlation_clustermap_task import (
    GeneticCorrelationClustermapTask,
    GeneticCorrSource,
    RGWithAsterix,
)

CT_LDSC_INITIAL_PLOT = GeneticCorrelationClustermapTask.create_std_with_clustering(
    "initial_genetic_correlation_by_ct_ldsc_plot",
    genetic_corr_source=GeneticCorrSource(task=CT_LDSC_INITIAL),
    plot_options=RGWithAsterix(color_scale="tropic"),
)
