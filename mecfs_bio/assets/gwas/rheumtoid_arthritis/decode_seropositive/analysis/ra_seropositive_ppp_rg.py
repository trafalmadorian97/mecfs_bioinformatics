from mecfs_bio.asset_generator.ukbb_ppp_rg_asset_generator import generate_ppp_rg_assets
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_standard_analysis import (
    SEROPOSITIVE_RA_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import PppRgConfig

SEROPOSITIVE_RA_PPP_RG_TASKS = generate_ppp_rg_assets(
    "seropositive_ra_ppp_rg",
    trait_task=SEROPOSITIVE_RA_STANDARD_ANALYSIS.assign_rsids_task_group.join_task,
    trait_pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
)

SEROPOSITIVE_RA_PPP_RG_TASKS_CIS_EXCLUDED = generate_ppp_rg_assets(
    "seropositive_ra_ppp_rg_cis_excluded",
    trait_task=SEROPOSITIVE_RA_STANDARD_ANALYSIS.assign_rsids_task_group.join_task,
    trait_pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
    config=PppRgConfig(
        variant_set="cis_excluded"
    ),
)
