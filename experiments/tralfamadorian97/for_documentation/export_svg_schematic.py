from mecfs_bio.assets.diagrams.coloc_configuration_grid import COLOC_CONFIGURATION_GRID_DIAGRAM
from mecfs_bio.assets.diagrams.ld_same_branch_correlation import LD_SAME_BRANCH_CORRELATION_DIAGRAM
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures([COLOC_CONFIGURATION_GRID_DIAGRAM])
    push_figures()

if __name__ == '__main__':
    go()