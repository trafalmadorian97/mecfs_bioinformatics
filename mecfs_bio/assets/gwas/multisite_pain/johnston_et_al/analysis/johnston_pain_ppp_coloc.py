"""
Exploratory colocalization sweep of Johnston et al.'s multisite chronic pain GWAS (quantitative)
with the UKB-PPP proteins. The trait standard deviation is estimated from allele frequency and N.
The lead-variant threshold is relaxed to 1e-5 to widen the candidate set.
"""

from mecfs_bio.asset_generator.ppp_coloc_asset_generator import generate_ppp_coloc_tasks
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.processed.johnston_pain_build_38_harmonized import (
    JOHNSTON_PAIN_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import QuantitativeTraitModel
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
)

JOHNSTON_PAIN_PPP_COLOC_FULL = generate_ppp_coloc_tasks(
    base_name="johnston_et_al_pain_ppp_coloc_full",
    trait_task=JOHNSTON_PAIN_BUILD_38_HARMONIZED,
    trait_model=QuantitativeTraitModel(sd_y=None),
    config=ColocRunConfig(
        exclude_mhc=True,
        selection=SelectionConfig(
            max_pairs=4000,
            max_per_locus=200,
        ),
        lead_p_threshold=1e-5,
    ),
)
