from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.polyfun.precomputed_prior.polyfun_precomputed_prior import \
    POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHT_CHR_8_22


def go():
    DEFAULT_RUNNER.run([

                        POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHT_CHR_8_22
                        ])

if __name__ == '__main__':
    go()