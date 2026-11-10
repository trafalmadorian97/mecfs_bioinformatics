"""
Driver to colocalize each trait with the Western et al. 2024 CSF aptamers: raw results and
display table per trait.

Alzheimer's runs first as the positive control (main() runs only it); the other traits and
the exploratory DecodeME runs are started explicitly once its results have been checked
against Western et al.'s reported AD proteins.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.protein_coloc_asset_generator import ColocTaskGroup
from mecfs_bio.assets.gwas.als.van_rheenen_et_al_2021.analysis.van_rheenen_2021_als_csf_coloc import (
    VAN_RHEENEN_2021_ALS_CSF_COLOC,
)
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.analysis.bellenguez_alzheimers_csf_coloc import (
    BELLENGUEZ_ALZHEIMERS_CSF_COLOC,
)
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.analysis.kerrebijn_fibro_csf_coloc import (
    KERREBIJN_FIBRO_CSF_COLOC,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.analysis.liu_et_al_2023_ibd_csf_coloc import (
    LIU_ET_AL_2023_IBD_CSF_COLOC,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_csf_coloc import (
    WILLER_LDL_CSF_COLOC,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_csf_coloc import (
    DECODE_ME_CSF_COLOC,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_csf_coloc_exploratory import (
    DECODE_ME_CSF_COLOC_LEAD_P_1E_5,
    DECODE_ME_CSF_COLOC_LEAD_P_1E_6,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_pain_csf_coloc import (
    JOHNSTON_PAIN_CSF_COLOC,
)
from mecfs_bio.assets.gwas.parkinsons.nalls_et_al_2019.analysis.nalls_2019_parkinsons_csf_coloc import (
    NALLS_2019_PARKINSONS_CSF_COLOC,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_csf_coloc import (
    RA_SEROPOSITIVE_CSF_COLOC,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.pgc_2022_sch_csf_coloc import (
    PGC_2022_SCH_CSF_COLOC,
)


def _run(group: ColocTaskGroup) -> None:
    DEFAULT_RUNNER.run(
        group.terminal_tasks(), incremental_save=True, must_rebuild_transitive=[]
    )


def run_bellenguez_alzheimers_csf_coloc() -> None:
    _run(BELLENGUEZ_ALZHEIMERS_CSF_COLOC)


def run_decode_me_csf_coloc() -> None:
    _run(DECODE_ME_CSF_COLOC)


def run_decode_me_csf_coloc_exploratory() -> None:
    _run(DECODE_ME_CSF_COLOC_LEAD_P_1E_6)
    _run(DECODE_ME_CSF_COLOC_LEAD_P_1E_5)


def run_ra_seropositive_csf_coloc() -> None:
    _run(RA_SEROPOSITIVE_CSF_COLOC)


def run_willer_ldl_csf_coloc() -> None:
    _run(WILLER_LDL_CSF_COLOC)


def run_liu_et_al_2023_ibd_csf_coloc() -> None:
    _run(LIU_ET_AL_2023_IBD_CSF_COLOC)


def run_pgc_2022_sch_csf_coloc() -> None:
    _run(PGC_2022_SCH_CSF_COLOC)


def run_nalls_2019_parkinsons_csf_coloc() -> None:
    _run(NALLS_2019_PARKINSONS_CSF_COLOC)


def run_van_rheenen_2021_als_csf_coloc() -> None:
    _run(VAN_RHEENEN_2021_ALS_CSF_COLOC)


def run_pain_csf_coloc() -> None:
    _run(KERREBIJN_FIBRO_CSF_COLOC)
    _run(JOHNSTON_PAIN_CSF_COLOC)


def main() -> None:
    run_bellenguez_alzheimers_csf_coloc()


if __name__ == "__main__":
    main()
