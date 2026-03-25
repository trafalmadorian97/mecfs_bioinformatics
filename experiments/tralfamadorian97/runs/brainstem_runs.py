from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.brainstem.whole_brainstem.xue_et_al.analysis.xue_whole_brainstem_standard_analysis import \
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.brainstem.whole_brainstem.xue_et_al.raw.raw_xue_whole_brainstem import \
    XUE_WHOLE_BRAINSTEM_VOLUME_RAW


def run_bs():
    DEFAULT_RUNNER.run(
[

    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.hba_magma_tasks.magma_hba_result_plot_task

]
            # XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.get_terminal_tasks()
            # XUE_WHOLE_BRAINSTEM_VOLUME_RAW,

         ,
        must_rebuild_transitive=[
            XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.hba_magma_tasks.magma_hba_result_plot_task
        ],
        incremental_save=True
    )


if __name__ == "__main__":
    run_bs()
