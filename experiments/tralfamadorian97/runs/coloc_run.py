from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_ppp_coloc import \
    RA_SEROPOSITIVE_COLOC_FULL


def go():
    DEFAULT_RUNNER.run([
        RA_SEROPOSITIVE_COLOC_FULL.display_task
    ])

if __name__ == '__main__':
    go()