"""
Colocalization of Nalls et al. 2019's Parkinson's disease GWAS (case-control, without 23andMe)
with the Western et al. 2024 CSF aptamers.

Case fraction 0.0698 = (15,056 + 18,618 proxy cases) / 482,730, counting the UK Biobank proxy
cases as cases (proxies dilute effects slightly).
"""

from mecfs_bio.asset_generator.csf_coloc_asset_generator import generate_csf_coloc_tasks
from mecfs_bio.assets.gwas.parkinsons.nalls_et_al_2019.processed.nalls_2019_build_38_harmonized import (
    NALLS_2019_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel

NALLS_2019_PARKINSONS_CSF_COLOC = generate_csf_coloc_tasks(
    base_name="nalls_2019_parkinsons",
    trait_task=NALLS_2019_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(case_fraction=0.0698),
)
