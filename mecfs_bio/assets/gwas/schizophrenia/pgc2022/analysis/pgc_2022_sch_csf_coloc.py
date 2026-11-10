"""Colocalization of the PGC 2022 schizophrenia GWAS (case-control) with the Western et al. 2024 CSF aptamers."""

from mecfs_bio.asset_generator.csf_coloc_asset_generator import generate_csf_coloc_tasks
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.auxiliary.prevalance_info import (
    PGC_2022_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.pgc_2022_sch_build_38_harmonized import (
    PGC_2022_SCH_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel

PGC_2022_SCH_CSF_COLOC = generate_csf_coloc_tasks(
    base_name="pgc_2022_sch",
    trait_task=PGC_2022_SCH_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(
        case_fraction=PGC_2022_PREVALENCE_INFO.sample_prevalence
    ),
)
