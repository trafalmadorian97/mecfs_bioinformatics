from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.partitioned_model_baseline_ld_scores_v1_1 import (
    THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES_V1_1,
)
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

BASE_MODEL_PARTITIONED_LD_SCORES_V1_1_EXTRACTED = ExtractTarGzipTask.create(
    asset_id="base_model_partitioned_ld_scores_v1_1_extracted",
    source_task=THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES_V1_1,
    read_mode="r",
)
