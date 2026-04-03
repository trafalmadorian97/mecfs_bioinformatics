from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.mv_ldl_heritability_task import MV_LDL_HERITABILITY_TASK, \
    MV_LDL_LDSC_RESULTS_MARKDOWN


def run_miscl_analysis():
    DEFAULT_RUNNER.run(
        [
            MV_LDL_LDSC_RESULTS_MARKDOWN
        # MV_LDL_HERITABILITY_TASK
            ],
        # [MILLION_VETERAN_MI_EUR_DATA_RAW]+MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.get_terminal_tasks(),
        # [MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task],
        # [MI_EUR_MANHATTAN],

        incremental_save=True,
        must_rebuild_transitive=[

            MV_LDL_LDSC_RESULTS_MARKDOWN
        ]
        # must_rebuild_transitive=[MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task]


    )


if __name__ == "__main__":
    run_miscl_analysis()
