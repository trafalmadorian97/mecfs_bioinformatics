"""Annotation coefficients of the external (precomputed PolyFun) prior against
the internal (DecodeME L2-regularized S-LDSC) prior."""

from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    PRECOMPUTED_POLYFUN_PRIOR_SOURCE,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.decode_me_l2_sldsc_prior import (
    DECODE_ME_L2_SLDSC_PRIOR_SOURCE,
)
from mecfs_bio.build_system.task.annotation_weights.annotation_weights_comparison_task import (
    AnnotationWeightsComparisonTask,
)

DECODE_ME_ANNOTATION_WEIGHTS_COMPARISON = AnnotationWeightsComparisonTask.create(
    asset_id="decode_me_polyfun_prior_annotation_weights_comparison",
    external_weights_tasks=PRECOMPUTED_POLYFUN_PRIOR_SOURCE.distinct_explanation_weights_tasks(),
    internal_weights_tasks=DECODE_ME_L2_SLDSC_PRIOR_SOURCE.distinct_explanation_weights_tasks(),
)
