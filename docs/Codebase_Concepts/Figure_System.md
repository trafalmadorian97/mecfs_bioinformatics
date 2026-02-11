# Figure System

Plots and tables are important to the presentation of data science results.  The `mecfs_bioinformatics` repo is transitioning to programmatic system for managing plots and tables (referred to in the codebase as "figures").  The components of this system are:

- The list `ALL_FIGURE_TASKS` in `mecfs_bio/figures/figure_tasks.py` specifies the tasks that should be used to generate figures.
- The script `mecfs_bio/figures/key_scripts/generate_figures.py` invokes the [build system](Build_System.md) to generate the assets correspond to the tasks in `ALL_FIGURE_TASKS`, then copies these figure assets into `docs/figs`.
- The script `mecfs_bio/figures/key_scripts/push_figures.py` uploads the figure assets in `docs/figs` to Github as a release.
- The script `mecfs_bio/figures/key_scripts/pull_figures.py` merges the contents of `docs/figs` with figures stored on Github.

## Standard Workflow


## Advantages