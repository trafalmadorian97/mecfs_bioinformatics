from mecfs_bio.assets.diagrams.ld_same_branch_correlation import LD_SAME_BRANCH_CORRELATION_DIAGRAM
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures([LD_SAME_BRANCH_CORRELATION_DIAGRAM])

if __name__ == '__main__':
    go()