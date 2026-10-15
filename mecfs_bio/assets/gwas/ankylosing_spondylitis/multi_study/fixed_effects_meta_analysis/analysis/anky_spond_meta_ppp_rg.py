from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.ankylosing_spondylitis.multi_study.fixed_effects_meta_analysis.processed.ank_spond_meta_analysis import (
    ANK_SPOND_FIXED_EFFECTS_META_ANALYSIS,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

ANK_SPOND_META_PPP_RG_CIS_EXCLUDED = generate_ppp_rg_assets(
    "ank_spond_ppp_rg_cis_excluded",
    trait_task=ANK_SPOND_FIXED_EFFECTS_META_ANALYSIS,
    # Trait N is read per variant from the meta-analysis's N column.
    config=PppRgConfig(variant_set="cis_excluded"),
)
