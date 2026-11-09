"""
Colocalization of Kerrebijn et al.'s European fibromyalgia GWAS (case-control) with the Western
et al. 2024 CSF aptamers. The sample sizes are effective sample sizes, so the case fraction is
0.5.
"""

from mecfs_bio.asset_generator.csf_coloc_asset_generator import generate_csf_coloc_tasks
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.processed.kerrebijn_fibro_build_38_harmonized import (
    KERREBIJN_FIBRO_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel

KERREBIJN_FIBRO_CSF_COLOC = generate_csf_coloc_tasks(
    base_name="kerrebijn_fibro",
    trait_task=KERREBIJN_FIBRO_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(case_fraction=0.5),
)
