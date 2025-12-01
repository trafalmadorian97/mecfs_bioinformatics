"""
Keep only the essential columns from the filtered genelist
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_filtered_gene_list_with_gene_metadata import (
    DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    CSVOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.log_markdown_pipe import LogMarkdownPipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA_DROP_COLS = PipeDataFrameTask.create(
    source_task=DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA,
    asset_id="decode_me_gwas_1_magma_filtered_gene_list_with_metadata_essential_cols",
    out_format=CSVOutFormat(sep=","),
    pipes=[
        SelectColPipe(
            cols_to_select=["GENE", "Gene name", "CHR", "P", "Gene description"],
        ),
        LogMarkdownPipe(),
    ],
)
