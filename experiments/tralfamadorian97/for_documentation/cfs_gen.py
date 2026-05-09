from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_ldsc import \
    DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD
from mecfs_bio.assets.gwas.me_cfs.million_veterans.analysis.million_veterans_cfs_standard_analysis import \
    MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures


def go():

    generate_figures([
        # MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP.heritability_markdown_task_unwrap
        DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD
    ])


if __name__ == '__main__':
    go()

