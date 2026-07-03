from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_standard_analysis import \
    DECODE_SEROPOSITIVE_RA_STANDARD_ANALYSIS, \
    SEROPOS_RA_FILTERED_FOR_FREQ
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.raw.download_raw_ra_seropositive import \
    DECODE_SEROPOSITIVE_RA_RAW


def go():
    DEFAULT_RUNNER.run(
        targets=DECODE_SEROPOSITIVE_RA_STANDARD_ANALYSIS.terminal_tasks(),
        # deCODE EAFrq was being read as a percentage (0-100) instead of a fraction.
        # The fix inserts an EAFrq /100 scaling step upstream of the frequency filter,
        # so the filter and everything downstream of it (sumstats, harmonization, MAGMA,
        # S-LDSC, cell-type LDSC) must be regenerated from corrected frequencies. Force a
        # transitive rebuild from the frequency filter; the new scaling task builds fresh
        # as its dependency, and the raw download / parquet conversion are reused.
        must_rebuild_transitive=[
            SEROPOS_RA_FILTERED_FOR_FREQ,
        ],
    )


if __name__ == '__main__':
    go()