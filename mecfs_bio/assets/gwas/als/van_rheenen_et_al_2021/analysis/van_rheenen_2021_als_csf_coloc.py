"""
Colocalization of van Rheenen et al. 2021's European ALS GWAS (case-control) with the Western
et al. 2024 CSF aptamers. The sample sizes are effective sample sizes, so the case fraction is
0.5.
"""

from mecfs_bio.asset_generator.csf_coloc_asset_generator import generate_csf_coloc_tasks
from mecfs_bio.assets.gwas.als.van_rheenen_et_al_2021.processed.van_rheenen_2021_build_38_harmonized import (
    VAN_RHEENEN_2021_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel

VAN_RHEENEN_2021_ALS_CSF_COLOC = generate_csf_coloc_tasks(
    base_name="van_rheenen_2021_als_eur",
    trait_task=VAN_RHEENEN_2021_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(case_fraction=0.5),
)
