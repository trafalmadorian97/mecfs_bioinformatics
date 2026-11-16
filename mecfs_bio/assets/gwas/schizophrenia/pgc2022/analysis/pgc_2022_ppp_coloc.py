from mecfs_bio.asset_generator.ppp_coloc_asset_generator import generate_ppp_coloc_tasks
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.auxiliary.prevalance_info import (
    PGC_2022_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.pgc_2022_sch_build_38_harmonized import (
    PGC_2022_SCH_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
)

PGC2022_PPP_COLOC_FULL = generate_ppp_coloc_tasks(
    base_name="decode_2022_sch_ppp_coloc_full",
    trait_task=PGC_2022_SCH_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(
        case_fraction=PGC_2022_PREVALENCE_INFO.sample_prevalence
    ),
    config=ColocRunConfig(
        exclude_mhc=True,
        selection=SelectionConfig(
            max_pairs=3000,
            max_per_locus=30,
        ),
        lead_p_threshold=1e-5,
    ),
)
