from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_combined_gene_lists import \
    DECODE_ME_GWAS_1_COMBINED_GENE_LISTS
from mecfs_bio.build_system.task.combine_gene_lists_task import ENSEMBL_ID_LABEL
from mecfs_bio.build_system.task.fetch_gget_info_task import FetchGGetInfoTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.sort_pipe import SortPipe

DECODE_ME_MASTER_GENE_LIST_WITH_GGET=FetchGGetInfoTask.create(
    asset_id="decode_me_gwas_1_master_gene_list_with_gget",
    source_df_task=DECODE_ME_GWAS_1_COMBINED_GENE_LISTS,
    ensembl_id_col=ENSEMBL_ID_LABEL,
    post_pipe=SortPipe([ENSEMBL_ID_LABEL]),
    out_format=ParquetOutFormat()
)