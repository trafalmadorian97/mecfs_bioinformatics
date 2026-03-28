from mecfs_bio.assets.gwas.multi_trait.polygenic_overlap.bivariate_mixer.mecfs_pain_bivariate_mixer import \
    MECFS_PAIN_BIVARIATE_MIXER
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures


def go():
    generate_figures([MECFS_PAIN_BIVARIATE_MIXER.result_table_markdown_task])


if __name__ == '__main__':
    go()

