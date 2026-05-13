from mecfs_bio.assets.gwas.height.yengo_2022.analysis.yengo_standard_analysis import YENGO_HEIGHT_STANDARD_ANALYSIS
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures


def go():
    generate_figures([
        YENGO_HEIGHT_STANDARD_ANALYSIS.heritability_markdown_task_unwrap

    ])


if __name__ == '__main__':
    go()

