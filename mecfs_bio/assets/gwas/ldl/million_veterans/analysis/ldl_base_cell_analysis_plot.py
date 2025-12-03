from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.ldl_base_cell_analysis_add_categories import (
    MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC_LABELED,
)
from mecfs_bio.build_system.task.sldsc_scatter_plot_task import SLDSCScatterPlotTask

LDL_BASE_CELL_PLOT = SLDSCScatterPlotTask.create(
    asset_id="ldl_base_cell_analysis_plot",
    source_task=MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC_LABELED,
)
