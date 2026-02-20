from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_counts_long import (
    YU_DRG_COUNTS_LONG,
)
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.raw.yu_drg_metadata_table import (
    YU_DRG_METADATA_TABLE,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.rename_col_by_position_pipe import (
    RenameColByPositionPipe,
)
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

YU_DRG_COUNTS_LONG_WITH_CELL_TYPE = JoinDataFramesTask.create_from_result_df(
    asset_id="yu_drg_counts_with_cell_types",
    result_df_task=YU_DRG_COUNTS_LONG,
    reference_df_task=YU_DRG_METADATA_TABLE,
    how="inner",
    left_on=["cell"],
    right_on=["cell"],
    out_format=ParquetOutFormat(),
    df_2_pipe=CompositePipe(
        [
            RenameColByPositionPipe(col_position=0, col_new_name="cell"),
            SelectColPipe(["cell", "cl.conserv_final"]),
        ]
    ),
)
