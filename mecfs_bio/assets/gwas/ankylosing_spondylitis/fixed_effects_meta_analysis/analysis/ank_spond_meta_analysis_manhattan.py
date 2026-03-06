from mecfs_bio.assets.gwas.ankylosing_spondylitis.fixed_effects_meta_analysis.processed.ank_spond_meta_analysis_sumstats import \
    ANK_SPOND_META_ANALYSIS_SUMSTATS
from mecfs_bio.build_system.task.gwaslab.gwaslab_manhattan_and_qq_plot_task import GWASLabManhattanAndQQPlotTask
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_mlog10p_pipe import ComputeMlog10pIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_p_from_beta_se import ComputePFromBetaSEPipeIfNeeded

ANK_SPOND_META_ANALYSIS_MANHATTAN=GWASLabManhattanAndQQPlotTask.create(
    ANK_SPOND_META_ANALYSIS_SUMSTATS,
    asset_id="ank_spond_meta_analysis_manhattan",
    plot_setting="m",
    pipe= CompositePipe([
ComputePFromBetaSEPipeIfNeeded(),
        ComputeMlog10pIfNeededPipe()
    ]),
)
