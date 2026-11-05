"""
Colocalization of deCODE's seropositive RA meta-GWAS (case-control) with the UKB-PPP proteins.
The case fraction is 0.5 because the sample sizes are effective sample sizes.
"""

from mecfs_bio.asset_generator.ppp_coloc_asset_generator import generate_ppp_coloc_tasks
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.processed.ra_seropositive_build_38_harmonized import (
    SEROPOS_RA_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel

RA_SEROPOSITIVE_PPP_COLOC = generate_ppp_coloc_tasks(
    base_name="decode_ra_seropositive",
    trait_task=SEROPOS_RA_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(case_fraction=0.5),
)
