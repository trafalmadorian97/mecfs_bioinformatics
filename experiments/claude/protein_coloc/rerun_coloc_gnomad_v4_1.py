"""
Reruns the four v1 trait x UKB-PPP colocalizations after switching genome-reference
harmonization (of the traits and of the protein windows) from the 1000 Genomes EUR panel to the
gnomAD v4.1 NFE panel, and after dropping the RA SNV-only filter.

Code changes do not invalidate the build cache, so the harmonized traits (and, for RA, the
gwaslab sumstats step whose filter changed) are forced to rebuild along with everything
downstream of them.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.analysis.liu_et_al_2023_ibd_ppp_coloc import (
    LIU_ET_AL_2023_IBD_PPP_COLOC,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_build_38_harmonized import (
    LIU_ET_AL_2023_IBD_EUR_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_ppp_coloc import (
    WILLER_LDL_PPP_COLOC,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.processed.willer_ldl_build_38_harmonized import (
    WILLER_LDL_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_ppp_coloc import (
    DECODE_ME_PPP_COLOC,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_38_harmonized import (
    DECODE_ME_GWAS_1_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_ppp_coloc import (
    RA_SEROPOSITIVE_PPP_COLOC,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.processed.ra_seropositive_build_38_harmonized import (
    SEROPOS_RA_BUILD_38_SUMSTATS,
)

_GROUPS_AND_FORCED = [
    (WILLER_LDL_PPP_COLOC, WILLER_LDL_BUILD_38_HARMONIZED),
    (DECODE_ME_PPP_COLOC, DECODE_ME_GWAS_1_BUILD_38_HARMONIZED),
    (RA_SEROPOSITIVE_PPP_COLOC, SEROPOS_RA_BUILD_38_SUMSTATS),
    (LIU_ET_AL_2023_IBD_PPP_COLOC, LIU_ET_AL_2023_IBD_EUR_BUILD_38_HARMONIZED),
]

if __name__ == "__main__":
    for group, forced in _GROUPS_AND_FORCED:
        DEFAULT_RUNNER.run(
            group.terminal_tasks(),
            incremental_save=True,
            must_rebuild_transitive=[forced],
        )
