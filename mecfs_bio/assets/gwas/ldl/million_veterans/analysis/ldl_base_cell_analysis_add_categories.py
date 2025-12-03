from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.ldl_base_cell_analysis_ldsc_filter_fdr import (
    MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC_FDR_FILTERED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.from_papers.finucane_2018_franke_gtex_categories import (
    FICUANE_2018_FRANKE_GTEX_CATEGORIES,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask

MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC_LABELED = JoinDataFramesTask.create_from_result_df(
    asset_id="million_veterans_ldl_base_cell_analysis_by_ldsc_fdr_filtered_categorized",
    result_df_task=MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC_FDR_FILTERED,
    reference_df_task=FICUANE_2018_FRANKE_GTEX_CATEGORIES,
    how="left",
    left_on=["Name"],
    right_on=["Tissue_Or_Cell"],
)
