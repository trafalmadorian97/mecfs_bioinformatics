"""
Markdown table of the DecodeME MAGMA gene-level results (with gene metadata),
for embedding in the docs in place of the hand-pasted table.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_filtered_gene_list_with_gene_metadata import (
    DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA,
)
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

DECODE_ME_GWAS_1_MAGMA_GENE_TABLE_MARKDOWN = (
    ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        source_task=DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA,
        asset_id="decode_me_gwas_1_magma_gene_table_markdown",
        pipe=SelectColPipe(
            cols_to_select=["GENE", "Gene name", "CHR", "P", "Gene description"],
        ),
    )
)
