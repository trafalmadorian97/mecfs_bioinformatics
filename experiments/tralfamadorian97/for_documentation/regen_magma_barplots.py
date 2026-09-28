from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_curated_gene_set_analysis import \
    DECODE_ME_CURATED_GENE_SET_ANALYSIS
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.magma_specific_tissue_bar_plot import \
    MAGMA_DECODE_ME_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [
            MAGMA_DECODE_ME_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
            DECODE_ME_CURATED_GENE_SET_ANALYSIS.bar_plot_task_full,

        ]
    )
    push_figures()

if __name__ == "__main__":
    go()