from mecfs_bio.asset_generator.ppp_coloc_asset_generator import generate_ppp_coloc_tasks
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.processed.bellenguez_build_38_harmonized import (
    BELLENGUEZ_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
)

BELLENGUEZ_PPP_COLOC_FULL = generate_ppp_coloc_tasks(
    base_name="bellenguez_alzheimers_ppp_coloc_full",
    trait_task=BELLENGUEZ_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(case_fraction=0.1763),
    config=ColocRunConfig(
        exclude_mhc=True,
        selection=SelectionConfig(
            max_pairs=4000,
            max_per_locus=30,
        ),
        lead_p_threshold=1e-5,
    ),
)
