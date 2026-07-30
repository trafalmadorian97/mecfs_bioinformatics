from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.migraine.multistudy.ct_ldsc_migraine_studies import MIGRAINE_CROSS_STUDY_CT_LDSC_ASSET_GENERATOR


def go():
    DEFAULT_RUNNER.run(
        MIGRAINE_CROSS_STUDY_CT_LDSC_ASSET_GENERATOR.terminal_tasks()
    )

if __name__ == '__main__':
    go()