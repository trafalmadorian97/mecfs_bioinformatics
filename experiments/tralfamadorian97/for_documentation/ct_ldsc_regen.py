from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_fast_plot import CT_LDSC_FAST_INITIAL_PLOT
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [
            CT_LDSC_FAST_INITIAL_PLOT
        ]
    )

if __name__ == '__main__':
    go()