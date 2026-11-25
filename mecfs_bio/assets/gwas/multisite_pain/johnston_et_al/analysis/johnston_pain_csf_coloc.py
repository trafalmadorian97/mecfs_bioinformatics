"""
Colocalization of Johnston et al.'s multisite chronic pain GWAS (quantitative) with the Western
et al. 2024 CSF aptamers. The trait standard deviation is estimated from allele frequency and N.
"""

from mecfs_bio.asset_generator.csf_coloc_asset_generator import generate_csf_coloc_tasks
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.processed.johnston_pain_build_38_harmonized import (
    JOHNSTON_PAIN_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import QuantitativeTraitModel

JOHNSTON_PAIN_CSF_COLOC = generate_csf_coloc_tasks(
    base_name="johnston_et_al_pain",
    trait_task=JOHNSTON_PAIN_BUILD_38_HARMONIZED,
    trait_model=QuantitativeTraitModel(sd_y=None),
)
