from mecfs_bio.assets.gwas.migraine.multistudy.ct_ldsc_migraine_studies import \
    MIGRAINE_CROSS_STUDY_CT_LDSC_ASSET_GENERATOR
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [MIGRAINE_CROSS_STUDY_CT_LDSC_ASSET_GENERATOR.aggregation_markdown_task]

    )

if __name__ == '__main__':
    go()