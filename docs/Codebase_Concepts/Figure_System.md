


# Figure System

Plots and tables are important to the presentation of data science results.  This section discusses how to add plots and tables (referred to as "figures") to documentation.  There are two main approaches:

- The link-based approach
- The programmatic approach



## Link-based approach

This is the simplest approach.

With this approach, if you wish to add a figure to a documentation page, simply save the figure as a `.png` file, upload it to a publicly accessible website (The repository wiki is convenient.  Ask a maintainer to be granted edit permission to the wiki), and then embed a reference in your documentation page:

```
![name_of_figure](https://<url_of_figure>.png)
```



## Programmatic approach

This approach is more complex, but more robust.

The components of the programmatic figure system are:

- The list `ALL_FIGURE_TASKS` in [figure_tasks.py][mecfs_bio.figures.figure_tasks] specifies the build-system [Tasks](Build_System.md#task) that generate programmatic figures.
- The script [generate_figures.py][mecfs_bio.figures.key_scripts.generate_figures] invokes the [build system](Build_System.md) to generate the assets corresponding to `ALL_FIGURE_TASKS`, then copies these figure assets into `docs/_figs`.
- The script [pull_figures.py][mecfs_bio.figures.key_scripts.pull_figures] downloads figures from Github and merges them with the contents of `docs/_figs`.
- The script [push_figures.py][mecfs_bio.figures.key_scripts.push_figures] uploads the figure assets in `docs/figs` to Github as a release. Running this script may require permission from a repository maintainer.


### Standard Workflow

Suppose that you have analyzed a genomic dataset and generated figures.  You wish to publish these figures and an associated write-up to the project documentation page.  Follow these steps:

- Add the Tasks that generate your figures to `ALL_FIGURE_TASKS`.
- Either invoke the [generate_figures.py][mecfs_bio.figures.key_scripts.generate_figures] script to generate all figures, or write your own script to call the [generate_figures][mecfs_bio.figures.key_scripts.generate_figures.generate_figures] function on just your newly added Tasks. In either case, your figures will be copied into `docs/_figs`
- Document your analysis by adding a markdown file to `docs/analysis`.  In your write-up, include your figures by referencing their location in `docs/_figs`.
- Upload your figures using [push_figures.py][mecfs_bio.figures.key_scripts.push_figures].
- Create a pull request with your changes (see [Standard Workflow](../Getting_Started/b_Standard_Workflow.md)).


### Advantages


This programmatic approach is undoubtedly more complex than the link-based approach. Nevertheless, this workflow, in which figures are represented as Tasks, has several advantages;

- **Figure Lineage**:  Since figures are represented as build-system Tasks, it is straightforward to determine exactly the datasets and analysis used to produce any given figure.  This supports that project principle of [automated reproducibility](../project_outline.md#principles) and is consistent with [data-science best practices](https://www.scientificdiscovery.dev/i/180950252/maximally-reproducible-charts).
- **Figure Upgrades**: Since figures are represented as build system Tasks, any improvement to that code that generates a particular class of figures can easily be propagated to all figures of that class.