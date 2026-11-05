"""Colocalization of the DecodeME GWAS-1 (case-control) with the UKB-PPP proteins."""

from mecfs_bio.asset_generator.ppp_coloc_asset_generator import generate_ppp_coloc_tasks
from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_38_harmonized import (
    DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel

DECODE_ME_PPP_COLOC = generate_ppp_coloc_tasks(
    base_name="decode_me_gwas_1",
    trait_task=DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(
        case_fraction=DECODE_ME_PREVALENCE_INFO.sample_prevalence
    ),
)
