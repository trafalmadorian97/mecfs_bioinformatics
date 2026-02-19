from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_extracted_counts_with_ensemlb_gene_names import (
    YU_DRG_COUNTS_WITH_ENSEMBL,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.unpivot_pipe import UnPivotPipe

YU_DRG_COUNTS_LONG = PipeDataFrameTask.create(
    YU_DRG_COUNTS_WITH_ENSEMBL,
    asset_id="yu_drg_counts_long",
    out_format=ParquetOutFormat(),
    pipes=[
        UnPivotPipe(
            on=None,
            index=["Gene stable ID"],
            variable_name="cell",
            value_name="expression",
        )
    ],
)
