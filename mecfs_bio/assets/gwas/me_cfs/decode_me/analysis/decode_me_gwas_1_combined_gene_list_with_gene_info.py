"""
Task to annotate the master gene list with gene symbols, names and UniProt protein information.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_combined_gene_lists import (
    DECODE_ME_GWAS_1_COMBINED_GENE_LISTS,
)
from mecfs_bio.assets.reference_data.gene_info.gene_info_reference import (
    GENE_INFO_REFERENCE,
)
from mecfs_bio.build_system.task.annotate_genes_from_reference_task import (
    AnnotateGenesFromReferenceTask,
)
from mecfs_bio.build_system.task.combine_gene_lists_task import ENSEMBL_ID_LABEL
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
)
from mecfs_bio.build_system.task.pipes.sort_pipe import SortPipe

DECODE_ME_MASTER_GENE_LIST_WITH_GENE_INFO = AnnotateGenesFromReferenceTask.create(
    asset_id="decode_me_gwas_1_master_gene_list_with_gene_info",
    source_df_task=DECODE_ME_GWAS_1_COMBINED_GENE_LISTS,
    ensembl_id_col=ENSEMBL_ID_LABEL,
    reference=GENE_INFO_REFERENCE,
    post_pipe=SortPipe([ENSEMBL_ID_LABEL]),
    out_format=ParquetOutFormat(),
)
