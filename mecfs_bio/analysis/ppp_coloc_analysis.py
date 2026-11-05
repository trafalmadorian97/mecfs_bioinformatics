"""
Driver to colocalize each v1 trait with the UKB-PPP proteins: raw results, display table and
stacked locus plots per trait.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.ppp_coloc_asset_generator import PppColocTaskGroup
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.analysis.liu_et_al_2023_ibd_ppp_coloc import (
    LIU_ET_AL_2023_IBD_PPP_COLOC,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_ppp_coloc import (
    WILLER_LDL_PPP_COLOC,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_ppp_coloc import (
    DECODE_ME_PPP_COLOC,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_ppp_coloc import (
    RA_SEROPOSITIVE_PPP_COLOC,
)


def _run(group: PppColocTaskGroup) -> None:
    DEFAULT_RUNNER.run(
        group.terminal_tasks(), incremental_save=True, must_rebuild_transitive=[]
    )


def run_willer_ldl_ppp_coloc() -> None:
    _run(WILLER_LDL_PPP_COLOC)


def run_decode_me_ppp_coloc() -> None:
    _run(DECODE_ME_PPP_COLOC)


def run_ra_seropositive_ppp_coloc() -> None:
    _run(RA_SEROPOSITIVE_PPP_COLOC)


def run_liu_et_al_2023_ibd_ppp_coloc() -> None:
    _run(LIU_ET_AL_2023_IBD_PPP_COLOC)


def main() -> None:
    run_willer_ldl_ppp_coloc()
    run_decode_me_ppp_coloc()
    run_ra_seropositive_ppp_coloc()
    run_liu_et_al_2023_ibd_ppp_coloc()


if __name__ == "__main__":
    main()
