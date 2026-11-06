"""
Exploratory: sub-genome-wide-significant DecodeME loci. Colocalization of the DecodeME GWAS-1
with the Western et al. 2024 CSF aptamers at relaxed lead p thresholds (1e-6 and 1e-5), to see
whether tentative DecodeME loci colocalize with CSF pQTLs. Results at these thresholds are
hypotheses, not findings.

DecodeME candidate loci (MHC excluded): 5 at p < 5e-8, 24 at p < 1e-6, 70 at p < 1e-5
(experiments/claude/csf_coloc/decode_me_relaxed_loci.py). The budget of 1,000 pairs, 20 per
locus, leaves about 14 pairs per locus even at 1e-5.
"""

from mecfs_bio.asset_generator.csf_coloc_asset_generator import generate_csf_coloc_tasks
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

_EXPLORATORY_SELECTION = SelectionConfig(max_pairs=1000, max_per_locus=20)
_TRAIT_MODEL = CaseControlTraitModel(
    case_fraction=DECODE_ME_PREVALENCE_INFO.sample_prevalence
)

DECODE_ME_CSF_COLOC_LEAD_P_1E_6 = generate_csf_coloc_tasks(
    base_name="decode_me_gwas_1_lead_p_1e_6",
    trait_task=DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
    trait_model=_TRAIT_MODEL,
    config=ColocRunConfig(lead_p_threshold=1e-6, selection=_EXPLORATORY_SELECTION),
)

DECODE_ME_CSF_COLOC_LEAD_P_1E_5 = generate_csf_coloc_tasks(
    base_name="decode_me_gwas_1_lead_p_1e_5",
    trait_task=DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
    trait_model=_TRAIT_MODEL,
    config=ColocRunConfig(lead_p_threshold=1e-5, selection=_EXPLORATORY_SELECTION),
)
