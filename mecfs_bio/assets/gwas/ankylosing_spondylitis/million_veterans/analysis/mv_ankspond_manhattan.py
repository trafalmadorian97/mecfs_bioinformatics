from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_sumstats_harmonized import \
    MILLION_VETERANS_ANK_SPOND_SUMSTATS_HARMONIZED
from mecfs_bio.build_system.task.gwaslab.gwaslab_manhattan_and_qq_plot_task import GWASLabManhattanAndQQPlotTask
from mecfs_bio.build_system.task.pipes.compute_mlog10p_pipe import ComputeMlog10pIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePPipe

MV_ANK_SPOND_MANHATTAN=GWASLabManhattanAndQQPlotTask.create(
    MILLION_VETERANS_ANK_SPOND_SUMSTATS_HARMONIZED,
    asset_id="mv_ank_spond_manhattan",
    plot_setting="m",
    pipe= ComputeMlog10pIfNeededPipe(),
)
