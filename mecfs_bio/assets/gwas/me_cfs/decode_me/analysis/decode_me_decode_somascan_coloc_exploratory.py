"""
Exploratory: sub-genome-wide-significant DecodeME loci. Colocalization of the DecodeME GWAS-1
with the deCODE SomaScan SMP assays at a relaxed lead p threshold (1e-5), to see whether
tentative DecodeME loci colocalize with deCODE pQTLs. Results at this threshold are hypotheses,
not findings.

Loci in the extended MHC are tested too: its long-range LD strains coloc.abf's
single-causal-variant assumption, so an MHC colocalization here is weaker evidence still.

The budgets are set so that they never bind, and require_no_budget_drop makes the run fail if
one does: every candidate pair is tested. The raw results are also written with a flag for
pairs whose protein lead lies near a highly pleiotropic SMP pQTL (ST36).
"""

from mecfs_bio.asset_generator.decode_somascan_coloc_asset_generator import (
    generate_decode_somascan_coloc_tasks,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_38_harmonized import (
    DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.reference_data.research_paper_supplementary_material.eldjarn_et_al_2023.processed.eldjarn_et_al_2023_st36_pleiotropic_pqtls import (
    ELDJARN_ET_AL_2023_ST36_PLEIOTROPIC_PQTLS,
)
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import CaseControlTraitModel
from mecfs_bio.build_system.task.protein_coloc.coloc_pleiotropy_flag_task import (
    ColocPleiotropyFlagTask,
)
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
)

DECODE_ME_DECODE_SOMASCAN_SMP_COLOC_LEAD_P_1E_5 = generate_decode_somascan_coloc_tasks(
    base_name="decode_me_gwas_1_lead_p_1e_5",
    trait_task=DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
    trait_model=CaseControlTraitModel(
        case_fraction=DECODE_ME_PREVALENCE_INFO.sample_prevalence
    ),
    config=ColocRunConfig(
        lead_p_threshold=1e-5,
        exclude_mhc=False,
        selection=SelectionConfig(max_pairs=1_000_000, max_per_locus=1_000_000),
        require_no_budget_drop=True,
    ),
)

DECODE_ME_DECODE_SOMASCAN_SMP_COLOC_LEAD_P_1E_5_PLEIOTROPY = (
    ColocPleiotropyFlagTask.create(
        asset_id="decode_me_gwas_1_lead_p_1e_5_decode_somascan_smp_coloc_pleiotropy",
        coloc_task=DECODE_ME_DECODE_SOMASCAN_SMP_COLOC_LEAD_P_1E_5.coloc_task,
        pleiotropic_task=ELDJARN_ET_AL_2023_ST36_PLEIOTROPIC_PQTLS,
    )
)
