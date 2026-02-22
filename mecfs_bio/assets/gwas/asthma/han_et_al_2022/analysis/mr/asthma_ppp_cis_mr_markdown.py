"""
Produce a markdown table of significant proteins from MR analysis of Asthma GWAS
"""

from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.mr.asthma_ppp_cis_mr_label_with_uniprot import (
    ASTHMA_PP_CIS_MR_UNIPROT_LABELED,
)
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.build_system.task.pipes.sort_pipe import SortPipe
from mecfs_bio.build_system.task.two_sample_mr_task import (
    TSM_OUTPUT_B_COL,
    TSM_OUTPUT_P_COL,
    TSM_OUTPUT_SE_COL,
)
from mecfs_bio.constants.sun_et_al_pqtl_constants import (
    SUN_ASSAY_TARGET,
    SUN_TARGET_UNIPROT,
)
from mecfs_bio.constants.uniprot_constants import (
    UNIPROT_DAT_FUNCTION_COL,
)

ASTHMA_HAN_PPP_CIS_MR_MARKDOWN = (
    ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        ASTHMA_PP_CIS_MR_UNIPROT_LABELED,
        asset_id="asthma_han_cis_mr_markdown",
        pipe=CompositePipe(
            [
                SelectColPipe(
                    [
                        SUN_TARGET_UNIPROT,
                        SUN_ASSAY_TARGET,
                        TSM_OUTPUT_B_COL,
                        TSM_OUTPUT_SE_COL,
                        TSM_OUTPUT_P_COL,
                        UNIPROT_DAT_FUNCTION_COL,
                    ],
                ),
                SortPipe(by=[TSM_OUTPUT_P_COL]),
            ]
        ),
    )
)
