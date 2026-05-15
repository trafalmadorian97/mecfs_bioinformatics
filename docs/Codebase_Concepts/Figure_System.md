


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

- The list `ALL_FIGURE_TASKS` in {{ api_link("figure_tasks.py", "mecfs_bio.figures.figure_task_list") }} specifies the build-system [Tasks](Build_System.md#task) that generate programmatic figures.
- The script {{ api_link("generate_figures.py", "mecfs_bio.figures.key_scripts.generate_figures") }} invokes the [build system](Build_System.md) to generate the assets corresponding to `ALL_FIGURE_TASKS`, then copies these figure assets into `docs/_figs`[^regen_note].
- The committed manifest `mecfs_bio/figures/figures_manifest.json` maps every live figure path (relative to `docs/_figs`) to the SHA-256 of its contents. The manifest is the source of truth for which figures are part of the project.
- The script {{ api_link("pull_figures.py", "mecfs_bio.figures.key_scripts.pull_figures") }} reads the manifest and downloads any missing or out-of-date blobs from the figures GitHub release into `docs/_figs`. Pass `prune=True` to also delete local files not listed in the manifest.
- The script {{ api_link("push_figures.py", "mecfs_bio.figures.key_scripts.push_figures") }} hashes the contents of `docs/_figs`, updates the manifest, and uploads any blobs not yet present on the release. Each blob is stored on the release as a content-addressed asset (asset name = SHA-256), so pushes never overwrite existing blobs and concurrent updates by different collaborators surface as ordinary git merge conflicts on the manifest. Pass `prune=True` to drop manifest entries whose files are no longer present locally. Running this script may require permission from a repository maintainer.


### Standard Workflow

Suppose that you have analyzed a genomic dataset and generated figures.  You wish to publish these figures and an associated write-up to the project documentation page.  Follow these steps:

- Add the Tasks that generate your figures to `ALL_FIGURE_TASKS`.
- Run {{ api_link("publish_figures.py", "mecfs_bio.figures.key_scripts.publish_figures") }} script to generate all figures, copy them to `docs/_figs` and push them as part of a github release.
- Document your analysis by adding a markdown file to `docs/analysis`.  In your write-up, include your figures by referencing their location in `docs/_figs`.
- Commit the updated `figures_manifest.json` along with your other changes.
- Create a pull request with your changes (see [Standard Workflow](../Getting_Started/b_Standard_Workflow.md)).

To remove a figure, delete the corresponding Task from `ALL_FIGURE_TASKS`, then use `publish_figures.py`

### Advantages


This programmatic approach is undoubtedly more complex than the link-based approach. Nevertheless, this workflow, in which figures are represented as Tasks, has several advantages;

- **Figure Lineage**:  Since figures are represented as build-system Tasks, it is straightforward to determine exactly the datasets and analysis used to produce any given figure.  This supports that project principle of [automated reproducibility](../project_outline.md#principles) and is consistent with [data-science best practices](https://www.scientificdiscovery.dev/i/180950252/maximally-reproducible-charts).
- **Figure Upgrades**: Since figures are represented as build system Tasks, any improvement to that code that generates a particular class of figures can easily be propagated to all figures of that class.



[^regen_note]:  The script {{ api_link("regenerate_figures.py", "mecfs_bio.figures.key_scripts.regenerate_figures") }} is similar to `generate_figures`, but forces the tasks which create the figures to be rerun.  This useful if, for example, one of your plotting Tasks has changed and you wish to propagate this change to all plots generated with this Task.