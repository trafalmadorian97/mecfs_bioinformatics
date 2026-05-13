from mecfs_bio.assets.gwas.me_cfs.multistudy.analysis.genetic_correlation.ct_ldsc.ct_ldsc_mecfs_studies import (
    CFS_CT_LDSC_ASSET_GENERATOR,
)
from mecfs_bio.build_system.task.genetic_correlation_clustermap_task import (
    GeneticCorrelationClustermapTask,
    GeneticCorrSource,
    RGWithAsterisk,
)

CT_LDSC_CFS_CORR_PLOT = GeneticCorrelationClustermapTask.create_std_with_clustering(
    "genetic_correlation_across_cfs_studies",
    genetic_corr_source=GeneticCorrSource(
        task=CFS_CT_LDSC_ASSET_GENERATOR.aggregation_task
    ),
    plot_options=RGWithAsterisk(color_scale="tropic"),
)
