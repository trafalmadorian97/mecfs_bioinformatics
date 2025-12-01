"""
Add Gene symbols, Entrez IDs, and Gene descriptions to the list of significant genes produced by MAGMA.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_filtered_gene_list import (
    DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST,
)
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask

DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST_WITH_GENE_METADATA = (
    JoinDataFramesTask.create_from_result_df(
        asset_id="decode_me_gwas_1_magma_filtered_gene_list_with_metadata",
        result_df_task=DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST,
        reference_df_task=GENE_THESAURUS,
        how="left",
        left_on=["GENE"],
        right_on=["Gene stable ID"],
    )
)
