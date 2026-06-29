"""
Tasks to rename the individual S-LDSC result tables, and drop extraneous columns
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import (
    DECODE_ME_S_LDSC,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc_export.decode_me_sldsc_roadmap_add_modality_column import (
    DECODE_ME_S_LDSC_ASSAY_SPECIFIC,
    DROP_HP_COLS_PIPE,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    CSVOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe

DROP_LOWER_CASE_NAME_PIPE = DropColPipe(cols_to_drop=["lower case name"])

DECODE_ME_SLDSC_TABLES_NO_EXTRA_LABELS = [
    PipeDataFrameTask.create(
        source_task=DECODE_ME_S_LDSC.partitioned_tasks[ref].cell_analysis_task,
        pipes=[],
        out_format=CSVOutFormat(","),
        asset_id=f"decode_me_s_ldsc_{ref}_results",
        backend="polars",
    )
    for ref in ["cahoy_cns", "gtex_brain", "corces_atac"]
]

DECODE_ME_SLDSC_TABLES_EXTRA_WITH_EXTRA_LABELS = [
    PipeDataFrameTask.create(
        source_task=DECODE_ME_S_LDSC.partitioned_tasks[ref].add_categories_task_unwrap,
        pipes=[DROP_HP_COLS_PIPE, DROP_LOWER_CASE_NAME_PIPE],
        out_format=CSVOutFormat(","),
        asset_id=f"decode_me_s_ldsc_{ref}_results",
        backend="polars",
    )
    for ref in ["immgen"]
] + [
    PipeDataFrameTask.create(
        source_task=DECODE_ME_S_LDSC.partitioned_tasks[ref].add_categories_task_unwrap,
        pipes=[DROP_HP_COLS_PIPE],
        out_format=CSVOutFormat(","),
        asset_id=f"decode_me_s_ldsc_{ref}_results",
        backend="polars",
    )
    for ref in ["multi_tissue_gene_expression"]
]


DECODE_ME_SLDSC_ALL_TABLES = (
    DECODE_ME_SLDSC_TABLES_NO_EXTRA_LABELS
    + DECODE_ME_SLDSC_TABLES_EXTRA_WITH_EXTRA_LABELS
    + DECODE_ME_S_LDSC_ASSAY_SPECIFIC
)
