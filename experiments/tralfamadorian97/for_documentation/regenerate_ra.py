from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seronegative.analysis.ra_seronegative_standard_analysis import \
    SERONEGATIVE_RA_STANDARD_ANALYSIS
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [SERONEGATIVE_RA_STANDARD_ANALYSIS.tasks.magma_gene_manhattan_plot_unwrap]
    )
    push_figures()

if __name__ == '__main__':
    go()