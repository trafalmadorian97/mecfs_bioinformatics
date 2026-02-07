from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="rna_seq_data",
        sub_group="gtex",
        sub_folder=PurePath("raw"),
        filename="GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz",
        extension=".gz",
        id=AssetId("gtex_v10_rna_seq_median_tissue_expression_raw"),
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t", skip_rows=2)),
    ),
    url="https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz",
    md5_hash="13dfdd0cb73c5899cd0b102aee5ee6c3",
)
