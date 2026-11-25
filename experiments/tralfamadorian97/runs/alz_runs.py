from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.processed.bellenguez_build_38_harmonized import \
    BELLENGUEZ_BUILD_38_SUMSTATS


def go():
    DEFAULT_RUNNER.run([BELLENGUEZ_BUILD_38_SUMSTATS])

if __name__ == '__main__':
    go()