from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multi_trait.polygenic_overlap.initial_mixer_run import MULTISITE_PAIN_DECODE_ME_INITIAL_MIXER
from mecfs_bio.assets.reference_data.mixer.processed.mixer_g1000_plink_eur_extracted import MIXER_G1000_PLINK_EXTRACTED
from mecfs_bio.assets.reference_data.mixer.raw.mixer_g1000_plink_eur_raw import MIXER_RAW_G1000_PLINK_DATA


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
            MULTISITE_PAIN_DECODE_ME_INITIAL_MIXER
        ],
        incremental_save=True,
        must_rebuild_transitive=[
        ]
    )


if __name__ == "__main__":
    run_mixer()
