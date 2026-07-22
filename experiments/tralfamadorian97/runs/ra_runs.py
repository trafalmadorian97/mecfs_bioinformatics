from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_ppp_rg import \
    SEROPOSITIVE_RA_PPP_RG_TASKS, SEROPOSITIVE_RA_PPP_RG_TASKS_CIS_EXCLUDED
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_standard_analysis import \
    SEROPOSITIVE_RA_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.processed.ra_seropositive_preprocess import \
    SEROPOS_RA_FILTERED_FOR_FREQ
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.raw.download_raw_ra_seropositive import \
    DECODE_SEROPOSITIVE_RA_RAW
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import \
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE


def go():
    DEFAULT_RUNNER.run(
        targets=
            # SEROPOSITIVE_RA_STANDARD_ANALYSIS.tasks.heritability_markdown_task_unwrap
            # SEROPOSITIVE_RA_PPP_RG_TASKS.rg_task
            SEROPOSITIVE_RA_PPP_RG_TASKS_CIS_EXCLUDED.get_terminal_tasks()
        ,
        must_rebuild_transitive=[
            # THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE
            # SEROPOS_RA_FILTERED_FOR_FREQ,
        ],
    )


if __name__ == '__main__':
    go()