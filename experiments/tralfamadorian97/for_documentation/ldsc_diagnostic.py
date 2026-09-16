from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER

from mecfs_bio.assets.gwas.height.yengo_2022.analysis.yengo_standard_analysis import (
    YENGO_HEIGHT_STANDARD_ANALYSIS,
)
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures
from mecfs_bio.figures.key_scripts.push_figures import push_figures


def go():
    DEFAULT_RUNNER.run(
        [
            YENGO_HEIGHT_STANDARD_ANALYSIS.ldsc_diagnostic_plot_task_unwrap
        ],
        must_rebuild_transitive=[

            YENGO_HEIGHT_STANDARD_ANALYSIS.ldsc_diagnostic_plot_task_unwrap
        ]
    )

def go_fig():
    generate_figures([YENGO_HEIGHT_STANDARD_ANALYSIS.ldsc_diagnostic_plot_task_unwrap])
    push_figures()

if __name__ == '__main__':
    go_fig()