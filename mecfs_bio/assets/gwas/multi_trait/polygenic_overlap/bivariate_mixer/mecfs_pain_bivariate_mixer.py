from mecfs_bio.asset_generator.mixer_asset_generator import (
    bivariate_mixer_asset_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.mixer.decode_me_univariate_mixer import (
    DECODE_ME_UNIVARIATE_MIXER,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.mixer.johnston_et_al_univariate_mixer import (
    JOHNSTON_ET_AL_UNIVARIATE_MIXER,
)

MECFS_PAIN_BIVARIATE_MIXER = bivariate_mixer_asset_generator(
    base_name="initial_bivariate_mixer",
    trait_1_tasks=DECODE_ME_UNIVARIATE_MIXER,
    trait_2_tasks=JOHNSTON_ET_AL_UNIVARIATE_MIXER,
    apply_extract_to_test=True,  # Avoid using too much memory
    extra_fit_args=[],
)
