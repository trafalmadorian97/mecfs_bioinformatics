"""Colocalization of Willer et al.'s LDL GWAS (quantitative) with the Western et al. 2024 CSF aptamers."""

from mecfs_bio.asset_generator.csf_coloc_asset_generator import generate_csf_coloc_tasks
from mecfs_bio.assets.gwas.ldl.willer_et_al.processed.willer_ldl_build_38_harmonized import (
    WILLER_LDL_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import QuantitativeTraitModel

WILLER_LDL_CSF_COLOC = generate_csf_coloc_tasks(
    base_name="willer_ldl",
    trait_task=WILLER_LDL_BUILD_38_HARMONIZED,
    trait_model=QuantitativeTraitModel(sd_y=None),
)
