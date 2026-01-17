from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.mr.asthma_ppp_cis_two_sample_mr_extract import (
    HAN_2022_ASTHMA_CIS_PPP_TSMR_EXTRACT,
)
from mecfs_bio.build_system.task.fdr_multiple_testing_table_task import (
    MultipleTestingTableTask,
)
from mecfs_bio.build_system.task.two_sample_mr_task import TSM_OUTPUT_P_COL

HAN_2022_ASTHMA_CIS_PPP_TSMR_MULTIPLE_TESTING = (
    MultipleTestingTableTask.create_from_result_table_task(
        alpha=0.01,
        method="bonferroni",
        asset_id="asthma_two_sample_mr_cis_ppp_multiple_testing",
        p_value_column=TSM_OUTPUT_P_COL,
        source_task=HAN_2022_ASTHMA_CIS_PPP_TSMR_EXTRACT,
        apply_filter=True,
    )
)
