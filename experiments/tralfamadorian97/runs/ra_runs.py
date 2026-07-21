from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_ppp_rg import \
    SEROPOSITIVE_RA_PPP_RG_TASKS
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_standard_analysis import \
    SEROPOSITIVE_RA_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.processed.ra_seropositive_preprocess import \
    SEROPOS_RA_FILTERED_FOR_FREQ
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.raw.download_raw_ra_seropositive import \
    DECODE_SEROPOSITIVE_RA_RAW


def go():
    DEFAULT_RUNNER.run(
        targets=[
            # SEROPOSITIVE_RA_STANDARD_ANALYSIS.tasks.heritability_markdown_task_unwrap
            SEROPOSITIVE_RA_PPP_RG_TASKS.rg_task
        ],
        must_rebuild_transitive=[
            # SEROPOS_RA_FILTERED_FOR_FREQ,
        ],
    )


if __name__ == '__main__':
    go()