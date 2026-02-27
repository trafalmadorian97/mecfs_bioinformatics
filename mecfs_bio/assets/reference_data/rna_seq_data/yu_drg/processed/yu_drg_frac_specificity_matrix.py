"""
Create the cell/gene specificity matrix for the YU DRG dataset using the traditional fractional approach.

Here the specificity of gene X for cell Y is log(1+ expression of X in Y)/(sum_Z log(1+ expression of X in Z))
"""

from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_counts_long_with_cell_type import (
    YU_DRG_COUNTS_LONG_WITH_CELL_TYPE,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import CSVOutFormat
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.move_col_to_front_pipe import MoveColToFrontPipe
from mecfs_bio.build_system.task.pipes.pivot_pipe import PivotPipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.pipes.shifted_log_pipe import ShiftedLogPipeInclude
from mecfs_bio.build_system.task.specificity_frac_task import (
    NORMALIZED_MEAN,
    PrepareSpecificityFraction,
)
from mecfs_bio.constants.magma_constants import MAGMA_GENE_COL

YU_DRG_FRAC_SPECIFICITY_MATRIX = PrepareSpecificityFraction.create(
    asset_id="yu_drg_frac_specificity_matrix",
    long_count_df_task=YU_DRG_COUNTS_LONG_WITH_CELL_TYPE,
    cell_type_col="cl.conserv_final",
    count_col="expression",
    gene_col="Gene stable ID",
    cell_col="cell",
    min_cells_per_type=5,
    out_format=CSVOutFormat(sep="\t"),
    pre_pipe=ShiftedLogPipeInclude(
        base=10, pseudocount=1, cols_to_include=["expression"]
    ),
    post_pipe=CompositePipe(
        [
            PivotPipe(
                on="cl.conserv_final", index="Gene stable ID", values=NORMALIZED_MEAN
            ),
            RenameColPipe(old_name="Gene stable ID", new_name=MAGMA_GENE_COL),
            MoveColToFrontPipe(target_col=MAGMA_GENE_COL),
        ]
    ),
)
