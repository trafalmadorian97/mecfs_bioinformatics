"""Colocalization of the DecodeME GWAS-1 (case-control) with the UKB-PPP proteins."""

from mecfs_bio.asset_generator.ppp_coloc_asset_generator import generate_ppp_coloc_tasks
from mecfs_bio.asset_generator.protein_coloc_asset_generator import (
    generate_coloc_locus_stack_plot_task,
    generate_coloc_protein_info_task,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_38_harmonized import (
    DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
)

DECODE_ME_PPP_COLOC = generate_ppp_coloc_tasks(
    base_name="decode_me_gwas_1",
    trait_task=DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(
        case_fraction=DECODE_ME_PREVALENCE_INFO.sample_prevalence
    ),
)

DECODE_ME_PPP_COLOC_EXPLORATORY = generate_ppp_coloc_tasks(
    base_name="decode_me_gwas_1_exploratory",
    trait_task=DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(
        case_fraction=DECODE_ME_PREVALENCE_INFO.sample_prevalence
    ),
    config=ColocRunConfig(
        lead_p_threshold=1e-5,
        exclude_mhc=False,
        selection=SelectionConfig(max_pairs=200, max_per_locus=20),
    ),
)


DECODE_ME_PPP_COLOC_EXPLORATORY_2 = generate_ppp_coloc_tasks(
    base_name="decode_me_gwas_1_exploratory_2",
    trait_task=DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(
        case_fraction=DECODE_ME_PREVALENCE_INFO.sample_prevalence
    ),
    config=ColocRunConfig(
        lead_p_threshold=1e-5,
        exclude_mhc=False,
        selection=SelectionConfig(max_pairs=2000, max_per_locus=30),
    ),
)

DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS = generate_coloc_locus_stack_plot_task(
    asset_id="decode_me_gwas_1_exploratory_2_ppp_coloc_stack_plots",
    group=DECODE_ME_PPP_COLOC_EXPLORATORY_2,
)

DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO = generate_coloc_protein_info_task(
    asset_id="decode_me_gwas_1_exploratory_2_ppp_coloc_protein_info",
    group=DECODE_ME_PPP_COLOC_EXPLORATORY_2,
)
