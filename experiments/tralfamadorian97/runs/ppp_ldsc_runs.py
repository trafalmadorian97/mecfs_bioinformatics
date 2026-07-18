from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_heritability import PPP_SAMPLE_SIZES, \
    HAPMAP_3_PPP_HERITABILITY


def go():
    DEFAULT_RUNNER.run(
        [
            PPP_SAMPLE_SIZES,
            HAPMAP_3_PPP_HERITABILITY
        ]
    )


if __name__ == "__main__":
    go()