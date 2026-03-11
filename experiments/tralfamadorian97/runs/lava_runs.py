from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.lava.initial_lava import BASIC_G100_LAVA_ANALYSIS


def run_lava():
    DEFAULT_RUNNER.run(
        [
            # LAVA_UKBB_LD_REF
            # LAVA_G100_EUR_LD_REF_EXTRACTED,
            # LAVA_G1000_EUR_REF_LD,
            # DEFAULT_LAVA_LOCUS_FILE
            BASIC_G100_LAVA_ANALYSIS
        ],
        incremental_save=True,
        must_rebuild_transitive=[
        ]
    )


if __name__ == "__main__":
    run_lava()
