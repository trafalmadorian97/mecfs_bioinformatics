"""External (precomputed PolyFun) vs internal (DecodeME L2-regularized S-LDSC)
prior comparison at the DecodeME chr20:47,653,230 locus."""

from mecfs_bio.asset_generator.polyfun_prior_comparison_asset_generator import (
    generate_assets_polyfun_prior_comparison,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr20_47_653_230_l2_sldsc_prior import (
    POLYFUN_EXPLAIN_CHR20_47_L2_SLDSC_PRIOR,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr20_47_653_230 import (
    POLYFUN_EXPLAIN_CHR20_47,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    SecondaryPositionFromSnpid,
)

POLYFUN_PRIOR_COMPARISON_CHR20_47 = generate_assets_polyfun_prior_comparison(
    base_name="decode_me_polyfun_prior_comparison_chr20_47_653_230",
    external_prior=POLYFUN_EXPLAIN_CHR20_47,
    internal_prior=POLYFUN_EXPLAIN_CHR20_47_L2_SLDSC_PRIOR,
    secondary_position=SecondaryPositionFromSnpid(build_label="hg38"),
)
