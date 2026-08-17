from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.sbayesrc_gfwm.to_cojo_ma import (
    DECODE_ME_BUILD_37_ANNOVAR_RSID_COJO_MA,
)
from mecfs_bio.assets.reference_data.sbayes_rc_gwfm_ref.imputed_13m_v1.stage_reference_data_imputed_13m_v1 import (
    SBAYESRC_GWFM_IMPUTED13M_V1,
)
from mecfs_bio.build_system.task.sbayesrc.gctb_fine_map_task import GctbFineMapTask
from mecfs_bio.build_system.task.sbayesrc.gctb_gwfm_constants import GCTB_VERSION
from mecfs_bio.build_system.task.sbayesrc.gctb_mcmc_options import GctbMcmcOptions
from mecfs_bio.constants.gh_constants import GH_CONTAINER_REGISTRY

# A deliberately tiny MCMC chain: this run exists only to exercise the full remote
# pipeline end-to-end
DUMMY_GCTB_RUN = GctbFineMapTask.create(
    id="dummy_decode_me_gwas_1_annovar_rsid_sbayesrc_gwfm",
    ma_task=DECODE_ME_BUILD_37_ANNOVAR_RSID_COJO_MA,
    reference_task=SBAYESRC_GWFM_IMPUTED13M_V1,
    image=f"{GH_CONTAINER_REGISTRY}/gctb:{GCTB_VERSION}",
    mcmc_options=GctbMcmcOptions(
        chain_length=100,
        burn_in=50,
        thin=10,
        out_freq=10,
        seed=42,
    ),

)
