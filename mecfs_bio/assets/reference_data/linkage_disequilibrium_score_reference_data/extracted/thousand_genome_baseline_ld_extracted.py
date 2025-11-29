from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.partitioned_model_baseline_ld_scores import (
    THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES,
)
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED = (
    ExtractTarGzipTask.create_from_reference_file_task(
        asset_id="base_model_partitioned_ld_scores_extracted",
        source_task=THOUSAND_GENOME_BASE_MODEL_PARTITIONED_LD_SCORES,
        read_mode="r",
    )
)
