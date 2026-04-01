from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import \
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOL


def run_miscl_analysis():
    DEFAULT_RUNNER.run(
        # [MILLION_VETERAN_MI_EUR_DATA_RAW]+MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.get_terminal_tasks(),
        # [MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task],
        [THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOL],
        # [MI_EUR_MANHATTAN],
        incremental_save=True,
        must_rebuild_transitive=[THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOL]
        # must_rebuild_transitive=[MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task]


    )


if __name__ == "__main__":
    run_miscl_analysis()
