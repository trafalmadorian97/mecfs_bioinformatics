"""
Colocalization of Bellenguez et al. 2022's Alzheimer's disease GWAS (case-control) with the
Western et al. 2024 CSF aptamers. The positive control of the CSF coloc work: Western et al.
report AD causal CSF proteins.

Case fraction 0.1763 = (39,106 clinical + 46,828 proxy cases) / 487,511, counting the UK
Biobank proxy cases as cases. Western et al.'s CSF cohorts (Knight ADRC, ADNI, ...) likely
overlap the AD GWAS samples; coloc.abf assumes independent samples and this is not corrected.
"""

from mecfs_bio.asset_generator.csf_coloc_asset_generator import generate_csf_coloc_tasks
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.processed.bellenguez_build_38_harmonized import (
    BELLENGUEZ_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel

BELLENGUEZ_ALZHEIMERS_CSF_COLOC = generate_csf_coloc_tasks(
    base_name="bellenguez_alzheimers",
    trait_task=BELLENGUEZ_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(case_fraction=0.1763),
)
