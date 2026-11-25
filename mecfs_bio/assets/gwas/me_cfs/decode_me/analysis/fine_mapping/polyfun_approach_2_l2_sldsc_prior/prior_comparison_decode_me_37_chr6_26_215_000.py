"""External (precomputed PolyFun) vs internal (DecodeME L2-regularized S-LDSC)
prior comparison at the DecodeME chr6:26,215,000 locus."""

from mecfs_bio.asset_generator.polyfun_prior_comparison_asset_generator import (
    generate_assets_polyfun_prior_comparison,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr6_26_215_000_l2_sldsc_prior import (
    POLYFUN_EXPLAIN_CHR6_26_L2_SLDSC_PRIOR,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr6_26_215_000 import (
    POLYFUN_EXPLAIN_CHR6_26,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    SecondaryPositionFromSnpid,
)

POLYFUN_PRIOR_COMPARISON_CHR6_26 = generate_assets_polyfun_prior_comparison(
    base_name="decode_me_polyfun_prior_comparison_chr6_26_215_000",
    external_prior=POLYFUN_EXPLAIN_CHR6_26,
    internal_prior=POLYFUN_EXPLAIN_CHR6_26_L2_SLDSC_PRIOR,
    secondary_position=SecondaryPositionFromSnpid(build_label="hg38"),
)
