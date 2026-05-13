from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.mixer.decode_me_univariate_mixer import DECODE_ME_UNIVARIATE_MIXER
from mecfs_bio.assets.gwas.multi_trait.polygenic_overlap.bivariate_mixer.mecfs_pain_bivariate_mixer import \
    MECFS_PAIN_BIVARIATE_MIXER
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.mixer.johnston_et_al_univariate_mixer import \
    JOHNSTON_ET_AL_UNIVARIATE_MIXER
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures


def go():
    generate_figures([MECFS_PAIN_BIVARIATE_MIXER.result_table_markdown_task,
                      DECODE_ME_UNIVARIATE_MIXER.result_markdown_table_task,
                      JOHNSTON_ET_AL_UNIVARIATE_MIXER.result_markdown_table_task
                      ])


if __name__ == '__main__':
    go()

