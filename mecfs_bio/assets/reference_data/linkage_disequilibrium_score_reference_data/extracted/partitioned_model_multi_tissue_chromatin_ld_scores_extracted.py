from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.partitioned_model_multi_tissue_chromatin_ld_scores import (
    PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_RAW,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED = (
    ExtractTarGzipTask.create(
        asset_id="partitioned_model_multi_tissue_chromatin_ld_scores_extracted",
        source_task=PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_RAW,
        read_mode="r",
    )
)
PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED_RD = DiscardDepsWrapper(
    PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED
)
