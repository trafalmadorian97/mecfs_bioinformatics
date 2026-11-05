"""Colocalization of the Liu et al. 2023 European IBD GWAS (case-control) with the UKB-PPP proteins."""

from mecfs_bio.asset_generator.ppp_coloc_asset_generator import generate_ppp_coloc_tasks
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.auxiliary.prevalence_nfo import (
    LIU_ET_AL_IBD_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_build_38_harmonized import (
    LIU_ET_AL_2023_IBD_EUR_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel

LIU_ET_AL_2023_IBD_PPP_COLOC = generate_ppp_coloc_tasks(
    base_name="liu_et_al_2023_ibd_eur",
    trait_task=LIU_ET_AL_2023_IBD_EUR_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(
        case_fraction=LIU_ET_AL_IBD_PREVALENCE_INFO.sample_prevalence
    ),
)
