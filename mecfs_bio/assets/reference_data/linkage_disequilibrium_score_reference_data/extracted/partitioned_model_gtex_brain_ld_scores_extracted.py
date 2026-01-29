from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.partitioned_model_gtex_brain_ld_scores import (
    PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_RAW,
)
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_EXTRACTED = ExtractTarGzipTask.create(
    asset_id="partitioned_model_gtex_brain_ld_scores_extracted",
    source_task=PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_RAW,
    read_mode="r",
)
