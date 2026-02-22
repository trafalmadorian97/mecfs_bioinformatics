from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_lead_variants import (
    DECODE_ME_GWAS_1_LEAD_VARIANTS,
)
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

DECODE_ME_GWAS_1_GWASLAB_FILTERED_GENE_LIST_WITH_GENE_METADATA = (
    JoinDataFramesTask.create_from_result_df(
        asset_id="decode_me_gwas_1_gwaslab_filtered_gene_list_with_metadata",
        result_df_task=DECODE_ME_GWAS_1_LEAD_VARIANTS,
        reference_df_task=GENE_THESAURUS,
        how="left",
        left_on=["GENE"],
        right_on=["Gene name"],
        df_1_pipe=SelectColPipe(["GENE"]),
        df_2_pipe=SelectColPipe(
            ["Gene name", "NCBI gene (formerly Entrezgene) ID", "Gene stable ID"]
        ),
    )
)
