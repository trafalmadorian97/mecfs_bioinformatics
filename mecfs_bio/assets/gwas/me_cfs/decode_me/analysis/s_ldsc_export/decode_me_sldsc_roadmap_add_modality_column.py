"""
Split roadmap into different tables.  One per modality.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import (
    DECODE_ME_S_LDSC,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    CSVOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue

HYPOTHESIS_TESTING_COLUMNS = ["_Corrected P Value_", "Reject Null"]

DROP_HP_COLS_PIPE = DropColPipe(
    cols_to_drop=HYPOTHESIS_TESTING_COLUMNS,
)

DECODE_ME_S_LDSC_ASSAY_SPECIFIC = [
    PipeDataFrameTask.create(
        source_task=DECODE_ME_S_LDSC.partitioned_tasks[
            "multi_tissue_chromatin"
        ].add_categories_task_unwrap,
        asset_id=f"decode_me_s_ldsc_result_multi_tissue_chromatin_{assay}_only",
        backend="polars",
        out_format=CSVOutFormat(sep=","),
        pipes=[
            FilterRowsByValue(target_column="Epigenetic_Assay", valid_values=[assay]),
            DROP_HP_COLS_PIPE,
            DropColPipe(cols_to_drop=["Cell"]),
        ],
    )
    for assay in ["DNase", "H3K27ac", "H3K4me3", "H3K4me1", "H3K9ac", "H3K36me3"]
]
