"""External (precomputed PolyFun) vs internal (DecodeME L2-regularized S-LDSC)
prior comparison at the DecodeME chr15:54,925,638 locus."""

from mecfs_bio.asset_generator.polyfun_prior_comparison_asset_generator import (
    generate_assets_polyfun_prior_comparison,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.susie_explain_decode_me_37_chr15_54_925_638_l2_sldsc_prior import (
    POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr15_54_925_638 import (
    POLYFUN_EXPLAIN_CHR15_54,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    SecondaryPositionFromSnpid,
)

POLYFUN_PRIOR_COMPARISON_CHR15_54 = generate_assets_polyfun_prior_comparison(
    base_name="decode_me_polyfun_prior_comparison_chr15_54_925_638",
    external_prior=POLYFUN_EXPLAIN_CHR15_54,
    internal_prior=POLYFUN_EXPLAIN_CHR15_54_L2_SLDSC_PRIOR,
    secondary_position=SecondaryPositionFromSnpid(build_label="hg38"),
)
