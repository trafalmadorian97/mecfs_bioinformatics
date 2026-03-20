from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.mixer.decode_me_univariate_mixer import DECODE_ME_UNIVARIATE_MIXER
from mecfs_bio.assets.gwas.multi_trait.polygenic_overlap.bivariate_mixer.pain_mecfs_bivariate_mixer import \
    MECFS_PAIN_BIVARIATE_MIXER
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.mixer.johnston_et_al_univariate_mixer import \
    JOHNSTON_ET_AL_UNIVARIATE_MIXER
from mecfs_bio.assets.reference_data.mixer.raw.mixer_g1000_plink_eur_raw import MIXER_RAW_G1000_PLINK_DATA


def run_mixer():
    DEFAULT_RUNNER.run(
        # [
        # ]+DECODE_ME_UNIVARIATE_MIXER.terminal_tasks()
        #  JOHNSTON_ET_AL_UNIVARIATE_MIXER.terminal_tasks(),
        MECFS_PAIN_BIVARIATE_MIXER.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            # DECODE_ME_UNIVARIATE_MIXER.combine_task
        ]
    )


if __name__ == "__main__":
    run_mixer()
