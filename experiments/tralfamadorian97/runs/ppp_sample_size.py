from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap3_ppp_heritability import PPP_SAMPLE_SIZES


def go():
    DEFAULT_RUNNER.run(
        [
            PPP_SAMPLE_SIZES
        ]
    )


if __name__ == "__main__":
    go()