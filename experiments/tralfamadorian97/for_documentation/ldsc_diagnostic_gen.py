from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_ldsc_diagnostic_plot import \
    DECODE_ME_GWAS_1_LDSC_DIAGNOSTIC_PLOT
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [DECODE_ME_GWAS_1_LDSC_DIAGNOSTIC_PLOT]
    )
    push_figures()

if __name__ == '__main__':
    go()