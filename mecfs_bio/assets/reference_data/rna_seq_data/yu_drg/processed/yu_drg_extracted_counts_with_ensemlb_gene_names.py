from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
)
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_extracted_counts_data import (
    YU_DRG_EXTRACTED_COUNTS,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe

INDEX_COL_NAME = "__index_level_0__"

YU_DRG_COUNTS_WITH_ENSEMBL = JoinDataFramesTask.create_from_result_df(
    asset_id="yu_drg_counts_with_ensembl_id",
    result_df_task=YU_DRG_EXTRACTED_COUNTS,
    reference_df_task=GENE_THESAURUS,
    how="inner",
    left_on=[INDEX_COL_NAME],
    right_on=["Gene name"],
    out_format=ParquetOutFormat(),
    df_2_pipe=CompositePipe(
        [
            UniquePipe(by=["Gene name"], keep="first", order_by=["Gene stable ID"]),
            SelectColPipe(["Gene name", "Gene stable ID"]),
        ]
    ),
    out_pipe=DropColPipe([INDEX_COL_NAME]),
)
