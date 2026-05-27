from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_ldsc import \
    DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_magma_gene_plot import DECODE_ME_MAGMA_GENE_PLOT
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_magma_gene_plot_with_window import \
    DECODE_ME_MAGMA_GENE_PLOT_WITH_WINDOW
from mecfs_bio.assets.gwas.me_cfs.million_veterans.analysis.million_veterans_cfs_standard_analysis import \
    MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP
from mecfs_bio.assets.gwas.me_cfs.multistudy.analysis.genetic_correlation.ct_ldsc.ct_ldsc_mecfs_studies import \
    CFS_CT_LDSC_ASSET_GENERATOR
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():

    # generate_figures([
    #     # MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP.heritability_markdown_task_unwrap
    #     # DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD
    #     # MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP.magma_tasks.inner.bar_plot_task
    #     CFS_CT_LDSC_ASSET_GENERATOR.aggregation_markdown_task
    # ])
    regenerate_figures([
        DECODE_ME_MAGMA_GENE_PLOT,
        DECODE_ME_MAGMA_GENE_PLOT_WITH_WINDOW,
    ])


if __name__ == '__main__':
    go()

