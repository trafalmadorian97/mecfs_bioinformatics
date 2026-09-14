from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_liftover_to_37 import \
    DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37


def go():
    DEFAULT_RUNNER.run([DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37])


if __name__ == '__main__':
    go()