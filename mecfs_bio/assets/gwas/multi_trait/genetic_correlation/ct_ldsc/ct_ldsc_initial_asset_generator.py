"""
Task generator to use cross-trait LD Score Regression to estimate genetic correlation between a several traits.
Each pair of traits is a Task.  This facilitates caching of results
"""

from mecfs_bio.asset_generator.genetic_correlation_asset_generator import (
    genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_initial_sources import (
    CT_LDSC_INITIAL_SOURCES,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)

CT_LDSC_INITIAL_ASSET_GENERATOR = genetic_corr_by_ct_ldsc_asset_generator(
    "initial_rg",
    sources=CT_LDSC_INITIAL_SOURCES,
    ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
    build="19",
)
