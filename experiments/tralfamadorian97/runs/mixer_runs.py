from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.mixer.univariate_mixer import DECODE_ME_UNIVARIATE_MIXER


def run_mixer():
    DEFAULT_RUNNER.run(
        [
            # LAVA_UKBB_LD_REF
            # LAVA_G100_EUR_LD_REF_EXTRACTED,
            # LAVA_G1000_EUR_REF_LD,
            # DEFAULT_LAVA_LOCUS_FILE
            # BASIC_G100_LAVA_ANALYSIS
            # MIXER_RAW_G1000_PLINK_DATA,
            # MIXER_G1000_PLINK_EXTRACTED
            # MULTISITE_PAIN_DECODE_ME_INITIAL_MIXER
        ]+DECODE_ME_UNIVARIATE_MIXER.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            DECODE_ME_UNIVARIATE_MIXER.combine_task
        ]
    )


if __name__ == "__main__":
    run_mixer()
