from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_plot import CT_LDSC_INITIAL_PLOT
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.key_scripts.push_figures import push_figures
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures([CT_LDSC_INITIAL_PLOT])
    push_figures(figure_tasks=ALL_FIGURE_TASKS)


if __name__ == '__main__':
    go()