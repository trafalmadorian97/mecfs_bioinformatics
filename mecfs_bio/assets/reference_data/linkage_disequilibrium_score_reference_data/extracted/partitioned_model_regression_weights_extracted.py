from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.partitioned_model_regression_weights import \
    HM3_LD_SCORE_REGRESSION_WEIGHTS
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED = ExtractTarGzipTask.create_from_reference_file_task(
    asset_id="s_ldsc_hapmap3_regression_weights_extracted",
    source_task= HM3_LD_SCORE_REGRESSION_WEIGHTS ,
    read_mode="r"
)
