"""
Task to take combined list of candidate genes generated from MAGMA and GWASLAB and write out to a markdown file
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_combined_gene_list_with_gene_info import (
    DECODE_ME_MASTER_GENE_LIST_WITH_GENE_INFO,
)
from mecfs_bio.build_system.task.combine_gene_lists_task import ENSEMBL_ID_LABEL
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.constants.gene_info_constants import (
    GENE_INFO_NAME_COL,
    GENE_INFO_SUBCELLULAR_LOCALISATION_COL,
    GENE_INFO_SYMBOL_COL,
    GENE_INFO_UNIPROT_DESCRIPTION_COL,
    GENE_INFO_UNIPROT_ID_COL,
)

DECODE_ME_MASTER_GENE_LIST_AS_MARKDOWN = (
    ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        DECODE_ME_MASTER_GENE_LIST_WITH_GENE_INFO,
        "decode_me_gwas_1_combined_gene_list_markdown",
        pipe=SelectColPipe(
            [
                ENSEMBL_ID_LABEL,
                "sources",
                GENE_INFO_SYMBOL_COL,
                GENE_INFO_NAME_COL,
                GENE_INFO_UNIPROT_ID_COL,
                GENE_INFO_SUBCELLULAR_LOCALISATION_COL,
                GENE_INFO_UNIPROT_DESCRIPTION_COL,
            ]
        ),
    )
)
