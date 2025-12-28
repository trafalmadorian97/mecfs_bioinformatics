from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_combined_gene_list_with_gget import (
    DECODE_ME_MASTER_GENE_LIST_WITH_GGET,
)
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

DECODE_ME_MASTER_GENE_LIST_AS_MARKDOWN = (
    ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        DECODE_ME_MASTER_GENE_LIST_WITH_GGET,
        "decode_me_gwas_1_combined_gene_list_markdown",
        pipe=SelectColPipe(
            [
                "Ensembl ID",
                "sources",
                "uniprot_id",
                "primary_gene_name",
                "subcellular_localisation",
                "ensembl_description",
                "uniprot_description",
                "ncbi_description",
            ]
        ),
    )
)
