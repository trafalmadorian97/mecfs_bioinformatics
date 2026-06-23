from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.build_38.decode_me_gwas_1_assign_rsids_build_38 import \
    DECODE_ME_GWAS_1_ASSIGN_RSIDS_BUILD_38


def go():
    DEFAULT_RUNNER.run(
        [DECODE_ME_GWAS_1_ASSIGN_RSIDS_BUILD_38],
        incremental_save=True,
    )


if __name__ == '__main__':
    go()