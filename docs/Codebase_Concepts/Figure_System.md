# Figure System

Plots and tables are important to the presentation of data science results.  The `mecfs_bioinformatics` repo is transitioning to programmatic system for managing plots and tables (referred to in the codebase as "figures").  The components of this system are:

- The list `ALL_FIGURE_TASKS` in [figure_tasks.py][mecfs_bio.figures.figure_tasks] specifies the build-system [Tasks](Build_System.md#task) that generate figures.
- The script [generate_figures.py][mecfs_bio.figures.key_scripts.generate_figures] invokes the [build system](Build_System.md) to generate the assets corresponding to `ALL_FIGURE_TASKS`, then copies these figure assets into `docs/_figs`.
- The script [pull_figures.py][mecfs_bio.figures.key_scripts.pull_figures] downloads figures from Github and merges them with the contents of `docs/_figs`.
- The script [push_figures.py][mecfs_bio.figures.key_scripts.push_figures] uploads the figure assets in `docs/figs` to Github as a release. Running this script may require permission from a repository maintainer.


## Standard Workflow

Suppose that you have run some analysis on a genomic dataset and produced some figures.  You wish to publish these figures and an associated write-up to the project documentation page.  Follow these steps:

- Add the tasks that generate your figures to `ALL_FIGURE_TASKS`.
- Either invoke the [generate_figures.py][mecfs_bio.figures.key_scripts.generate_figures] script to generate all figures, or write your own script to call the [generate_figures][mecfs_bio.figures.key_scripts.generate_figures.generate_figures] function on just your newly added tasks. In either case, your figures will be copied into `docs/_figs`
- Document your analysis by adding a markdown file to `docs/analysis`.  In your write-up, include your figures by referencing their location in `docs/_figs`.
- Upload your figures using [push_figures.py][mecfs_bio.figures.key_scripts.push_figures].
- Create a pull request with your changes (see [Standard Workflow](../Getting_Started/b_Standard_Workflow.md)).


## Advantages


This programmatic workflow, in which figures are represented as Tasks has several advantages;

- **Figure Lineage**:  Since figures are represented as build-system tasks, it is straightforward to determine exactly the datasets and analysis used to produce any given figure.  This supports that project principle of [automated reproducibility](../project_outline.md#principles) and is consistent with [data-science best practices](https://www.scientificdiscovery.dev/i/180950252/maximally-reproducible-charts).
- **Figure Upgrades**: Since figures are represented as build system tasks, any improvement to that code that generates a particular class of figures can easily be propagated to all figures of that class.