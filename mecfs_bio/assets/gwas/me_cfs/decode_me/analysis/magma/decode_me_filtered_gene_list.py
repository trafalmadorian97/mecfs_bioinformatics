from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_ensembl_gene_analysis import (
    DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
)
from mecfs_bio.build_system.task.multiple_testing_table_task import (
    MultipleTestingTableTask,
)

DECODE_ME_GWAS_1_MAGMA_FILTERED_GENE_LIST = (
    MultipleTestingTableTask.create_from_magma_gene_analysis_task(
        asset_id="decode_me_gwas_1_magma_filtered_gene_list",
        alpha=0.01,
        method="fdr_bh",
        source_task=DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS,
    )
)
