from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.partitioned_model_baseline_ld_scores_v1_0 import (
    THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES_V1_0,
)
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

BASE_MODEL_PARTITIONED_LD_SCORES_V1_0_EXTRACTED = ExtractTarGzipTask.create(
    asset_id="base_model_partitioned_ld_scores_v1_0_extracted",
    source_task=THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES_V1_0,
    read_mode="r",
)
