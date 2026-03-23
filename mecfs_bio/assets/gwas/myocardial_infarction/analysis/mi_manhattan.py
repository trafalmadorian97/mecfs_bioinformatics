from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_standard_analysis import \
    MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS
from mecfs_bio.build_system.task.gwaslab.gwaslab_manhattan_and_qq_plot_task import GWASLabManhattanAndQQPlotTask
from mecfs_bio.build_system.task.pipes.compute_mlog10p_pipe import ComputeMlog10pIfNeededPipe

MI_EUR_MANHATTAN=GWASLabManhattanAndQQPlotTask.create(
    sumstats_task=MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task,
    asset_id="million_veterans_eur_mi_manhattan",
    plot_setting="m",
    pipe=ComputeMlog10pIfNeededPipe()
)