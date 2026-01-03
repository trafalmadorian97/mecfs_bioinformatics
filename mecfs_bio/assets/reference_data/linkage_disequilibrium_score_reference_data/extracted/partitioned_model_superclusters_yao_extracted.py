from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.partitioned_model_superclusters_yao import (
    PARTITIONED_MODEL_SUPERCLUSTERS_YAO,
)
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

PARTITIONED_MODEL_YAO_SUPERCLUSTERS_EXTRACTED = (
    ExtractTarGzipTask.create_from_reference_file_task(
        asset_id="s_ldsc_yao_superclusers_ld_scores_extracted",
        source_task=PARTITIONED_MODEL_SUPERCLUSTERS_YAO,
        read_mode="r",
    )
)
