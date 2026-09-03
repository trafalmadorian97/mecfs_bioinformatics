"""Reference asset: ridge surrogate weights (gamma_c) for the polyfun prior."""

from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotations import (
    BASELINE_LF_ANNOTATION_MATRIX,
)
from mecfs_bio.assets.reference_data.polyfun.precomputed_prior.polyfun_precomputed_prior import (
    COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    RidgeAnnotationWeightsTask,
)

BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS = RidgeAnnotationWeightsTask.create(
    asset_id="baseline_lf_2.2_ukb_annotation_ridge_weights",
    annotation_parquet_task=BASELINE_LF_ANNOTATION_MATRIX,
    snpvar_meta_task=COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
)
