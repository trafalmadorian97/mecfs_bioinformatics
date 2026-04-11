from mecfs_bio.assets.gwas.multi_trait.lcv.mi_lcv_analysis import MI_LCV_TASK_GROUP
from mecfs_bio.build_system.task.lcv.lcv_clustermap import LCVClustermapTask, LCVSource, GCPWithAsterisk

MI_LCP_PLOT=LCVClustermapTask.create_std_with_clustering(
    asset_id="mi_lcv_analysis_plot",
    source=LCVSource(
        task=MI_LCV_TASK_GROUP.agg_task,
    ),
    plot_options=GCPWithAsterisk()
)