from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.triglycerides.willer_et_al import \
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS
from mecfs_bio.build_system.scheduler.topological_scheduler import TopologicalSchedulerSettings


def run_miscl_analysis():
    DEFAULT_RUNNER.run(

            # WILLER_LDL_EUR_DATA_RAW
            # MV_LDL_LDSC_RESULTS_MARKDOWN
        # MV_LDL_HERITABILITY_TASK
        #     MI_LDL_WILLER_CORRELATION.terminal_tasks()
    # [
        # LDL_MI_LCV_ANALYSIS,
        WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS.get_terminal_tasks()
        # WILLER_TG_EUR_DATA_RAW
     # MECFS_PAIN_LCV_ANALYSIS
     # ]
            ,
        # +WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.get_terminal_tasks(),
        # [MILLION_VETERAN_MI_EUR_DATA_RAW]+MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.get_terminal_tasks(),
        # [MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task],
        # [MI_EUR_MANHATTAN],

        incremental_save=True,
        must_rebuild_transitive=[

            # MV_LDL_LDSC_RESULTS_MARKDOWN
        ],
        settings=TopologicalSchedulerSettings(
            print_progress=False
        )
        # must_rebuild_transitive=[MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task]


    )


if __name__ == "__main__":
    run_miscl_analysis()
