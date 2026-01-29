from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.partitioned_model_multi_tissue_gene_expr_ld_scores import (
    THOUSAND_GENOME_PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_RAW,
)
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED = ExtractTarGzipTask.create(
    asset_id="multi_tissue_gene_expression_partitioned_ld_scores_extracted",
    source_task=THOUSAND_GENOME_PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_RAW,
    read_mode="r",
)
