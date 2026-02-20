from mecfs_bio.assets.reference_data.rna_seq_data.gtex.gtex_v10_median_tissue_expression_rna_seq_raw import (
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ,
)
from mecfs_bio.build_system.task.extract_gzip_task import ExtractGzipTextFileTask

GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_EXTRACTED = (
    ExtractGzipTextFileTask.create_for_reference_file(
        source_file_task=GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ,
        asset_id="gtex_v10_rna_seq_median_tissue_expression_extracted",
    )
)
