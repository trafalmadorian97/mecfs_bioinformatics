from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.build_38.build_38_s_ldsc_export.build_38_roadmap_add_modaility import \
    BUILD_38_DECODE_ME_S_LDSC_ASSAY_SPECIFIC
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.build_38.decode_me_build_38_s_ldsc import DECODE_ME_S_LDSC_BUILD_38
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc_export.decode_me_sldsc_roadmap_add_modality_column import \
    DROP_HP_COLS_PIPE
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask, CSVOutFormat

BUILD_38_DECODE_ME_SLDSC_TABLES_NO_EXTRA_LABELS = [
    PipeDataFrameTask.create(
        source_task=DECODE_ME_S_LDSC_BUILD_38,
        pipes=[],
        out_format=CSVOutFormat(","),
        asset_id=f"build_38_decode_me_s_ldsc_{ref}_results",
        backend="polars",
    )
    for ref in ["cahoy_cns", "gtex_brain", "corces_atac"]
]

BUILD_38_DECODE_ME_SLDSC_TABLES_EXTRA_WITH_EXTRA_LABELS = [
    PipeDataFrameTask.create(
        source_task=DECODE_ME_S_LDSC_BUILD_38.partitioned_tasks[ref].add_categories_task_unwrap,
        pipes=[DROP_HP_COLS_PIPE],
        out_format=CSVOutFormat(","),
        asset_id=f"build_38_decode_me_s_ldsc_{ref}_results",
        backend="polars",
    )
    for ref in ["immgen", "multi_tissue_gene_expression"]
]


BUILD_38_DECODE_ME_SLDSC_ALL_TABLES = (
        BUILD_38_DECODE_ME_SLDSC_TABLES_NO_EXTRA_LABELS
        + BUILD_38_DECODE_ME_SLDSC_TABLES_EXTRA_WITH_EXTRA_LABELS
        + BUILD_38_DECODE_ME_S_LDSC_ASSAY_SPECIFIC
)
