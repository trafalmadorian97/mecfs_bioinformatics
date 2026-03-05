from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.processed.finngen_ank_spond_sumstats_harmonized import \
    FINGNEN_ANK_SPOND_SUMSTATS_HARMONIZED
from mecfs_bio.build_system.task.gwaslab.gwaslab_manhattan_and_qq_plot_task import GWASLabManhattanAndQQPlotTask

FINGEN_ANK_SPOND_MANHATTAN=GWASLabManhattanAndQQPlotTask.create(
FINGNEN_ANK_SPOND_SUMSTATS_HARMONIZED,
    asset_id="fingen_ank_spond_manhattan",
    plot_setting="m"
)