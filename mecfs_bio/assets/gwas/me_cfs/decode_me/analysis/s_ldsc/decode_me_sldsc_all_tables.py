from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc.decode_me_sldsc import DECODE_ME_S_LDSC
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.s_ldsc.decode_me_sldsc_roadmap_add_modality_column import \
    DECODE_ME_S_LDSC_ASSAY_SPECIFIC, DROP_HP_COLS_PIPE
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask, CSVOutFormat

DECODE_ME_SLDSC_TABLES_NO_EXTRA_LABELS=[
    DECODE_ME_S_LDSC.partitioned_tasks[ref].cell_analysis_task for ref in [
    "cahoy_cns",
    "gtex_brain",
    "corces_atac"
    ]
]

DECODE_ME_SLDSC_TABLES_EXTRA_WITH_EXTRA_LABELS=[
    PipeDataFrameTask.create(source_task=DECODE_ME_S_LDSC.partitioned_tasks[ref].add_categories_task_unwrap,
                      pipes=[
                          DROP_HP_COLS_PIPE
                      ],
                             out_format=CSVOutFormat(","),
                             asset_id=f"decode_me_s_ldsc_{ref}_results"
                      ) for ref in [
    "immgen",
    "multi_tissue_gene_expression"
]
]


DECODE_ME_SLDSC_ALL_TABLES = DECODE_ME_SLDSC_TABLES_NO_EXTRA_LABELS + DECODE_ME_SLDSC_TABLES_EXTRA_WITH_EXTRA_LABELS +DECODE_ME_S_LDSC_ASSAY_SPECIFIC