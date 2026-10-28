"""
Fast all-pairs CT-LDSC genetic correlation for the initial study traits, alongside the gwaslab
generator in ct_ldsc_initial_asset_generator for comparison.
"""

from mecfs_bio.asset_generator.fast_genetic_correlation_asset_generator import (
    fast_genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_initial_sources import (
    CT_LDSC_INITIAL_SOURCES,
)

CT_LDSC_FAST_INITIAL_ASSET_GENERATOR = fast_genetic_corr_by_ct_ldsc_asset_generator(
    "initial_rg_fast", sources=CT_LDSC_INITIAL_SOURCES
)
