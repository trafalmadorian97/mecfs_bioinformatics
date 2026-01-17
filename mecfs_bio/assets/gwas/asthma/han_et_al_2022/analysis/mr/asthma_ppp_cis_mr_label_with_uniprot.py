"""
Label significant proteins from MR analysis for Asthma GWAS using descriptions from UniProt
"""

from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.mr.asthma_ppp_cis_mr_multiple_testing import (
    HAN_2022_ASTHMA_CIS_PPP_TSMR_MULTIPLE_TESTING,
)
from mecfs_bio.assets.reference_data.uniprot.uniprot_lookup_table import UNIPROT_LOOKUP
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.str_split_exact_col import SplitExactColPipe
from mecfs_bio.build_system.task.two_sample_mr_task import TSM_OUTPUT_EXPOSURE_COL
from mecfs_bio.constants.sun_et_al_pqtl_constants import (
    SUN_ASSAY_TARGET,
    SUN_TARGET_UNIPROT,
)
from mecfs_bio.constants.uniprot_constants import UNIPROT_DAT_ID_COL

ASTHMA_PP_CIS_MR_UNIPROT_LABELED = JoinDataFramesTask.create_from_result_df(
    asset_id="asthma_ppp_cis_mr_labeled_with_uniprot",
    result_df_task=HAN_2022_ASTHMA_CIS_PPP_TSMR_MULTIPLE_TESTING,
    reference_df_task=UNIPROT_LOOKUP,
    how="left",
    df_1_pipe=SplitExactColPipe(
        TSM_OUTPUT_EXPOSURE_COL,
        split_by="_",
        new_col_names=(SUN_ASSAY_TARGET, SUN_TARGET_UNIPROT),
    ),
    left_on=[SUN_TARGET_UNIPROT],
    right_on=[UNIPROT_DAT_ID_COL],
    out_format=ParquetOutFormat(),
)
