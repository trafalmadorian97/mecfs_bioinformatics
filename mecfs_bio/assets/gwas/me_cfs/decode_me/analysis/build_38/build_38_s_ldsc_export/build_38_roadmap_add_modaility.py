from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.build_38.decode_me_build_38_s_ldsc import (
    DECODE_ME_S_LDSC_BUILD_38,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc_export.decode_me_sldsc_roadmap_add_modality_column import (
    DROP_HP_COLS_PIPE,
)
from mecfs_bio.build_system.task.dataframe_output import (
    CSVOutFormat,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue

BUILD_38_DECODE_ME_S_LDSC_ASSAY_SPECIFIC = [
    PipeDataFrameTask.create(
        source_task=DECODE_ME_S_LDSC_BUILD_38.partitioned_tasks[
            "multi_tissue_chromatin"
        ].add_categories_task_unwrap,
        asset_id=f"build_38_decode_me_s_ldsc_result_multi_tissue_chromatin_{assay}_only",
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
