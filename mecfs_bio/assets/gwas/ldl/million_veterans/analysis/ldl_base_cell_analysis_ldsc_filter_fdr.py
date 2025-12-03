from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.ldl_base_cell_analysis_by_ldsc import (
    MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC,
)
from mecfs_bio.build_system.task.fdr_multiple_testing_table_task import (
    FilterMultipleTestingTableTask,
)

MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC_FDR_FILTERED = (
    FilterMultipleTestingTableTask.create_from_result_table_task(
        alpha=0.05,
        p_value_column="Coefficient_P_value",
        method="fdr_bh",
        source_task=MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC,
        asset_id="million_veterans_ldl_base_cell_analysis_by_ldsc_fdr_filtered",
        apply_filter=False,
    )
)
