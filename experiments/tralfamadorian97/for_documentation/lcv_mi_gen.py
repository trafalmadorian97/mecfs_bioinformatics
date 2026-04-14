from mecfs_bio.assets.gwas.multi_trait.lcv.mi_lcv_analysis import MI_LCV_TASK_GROUP
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures


def go():
    generate_figures([
        MI_LCV_TASK_GROUP.downstream_trait_tables["MI"]
    ])


if __name__ == '__main__':
    go()

