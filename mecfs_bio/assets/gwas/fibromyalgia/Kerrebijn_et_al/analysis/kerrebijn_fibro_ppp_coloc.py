"""
Exploratory colocalization sweep of Kerrebijn et al.'s European fibromyalgia GWAS (case-control)
with the UKB-PPP proteins. The sample sizes are effective sample sizes, so the case fraction is
0.5. The lead-variant threshold is relaxed to 1e-5 to widen the candidate set.
"""

from mecfs_bio.asset_generator.ppp_coloc_asset_generator import generate_ppp_coloc_tasks
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.processed.kerrebijn_fibro_build_38_harmonized import (
    KERREBIJN_FIBRO_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
)

KERREBIJN_FIBRO_PPP_COLOC_FULL = generate_ppp_coloc_tasks(
    base_name="kerrebijn_fibro_ppp_coloc_full",
    trait_task=KERREBIJN_FIBRO_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(case_fraction=0.5),
    config=ColocRunConfig(
        exclude_mhc=True,
        selection=SelectionConfig(
            max_pairs=4000,
            max_per_locus=30,
        ),
        lead_p_threshold=1e-5,
    ),
)
