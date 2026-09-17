from mecfs_bio.assets.diagrams.liability_threshold_model import LIABILITY_THRESHOLD_MODEL_DIAGRAM
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures([
        LIABILITY_THRESHOLD_MODEL_DIAGRAM
    ])

if __name__ == '__main__':
    go()