from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_curated_gene_set_analysis import \
    DECODE_ME_CURATED_GENE_SET_ANALYSIS
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures


def go():
    generate_figures([
        DECODE_ME_CURATED_GENE_SET_ANALYSIS.bar_plot_task_full,
    ])


if __name__ == '__main__':
    go()

