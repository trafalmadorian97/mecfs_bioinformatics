# Build System

## Motivation
In data science, analysis pipelines often consist of many steps.  These steps can be a mixture of pure data-cleaning operations and more complex statistical modeling.  When the dataset under study is nontrivial in size, these steps can be quite slow.  This creates two challenges:

- **Iteration**: It is rare to run a pipeline once, produce an analysis, and be done. Usually, one must repeatedly tweak the steps, re-run the pipeline, and reexamine the result.  To avoid wasting time, it is therefore desirable that after each change, only impacted steps should be rerun.
- **Lineage**:  Given the complexity of many data science workflows, there is considerable room for error.  It is therefore desirable to be able to interrogate the final product of a workflow to trace its "lineage": the precise sequence of steps that produced it.

These challenges motivate the development of a data science build system.

## Framework 

The build system used in this project is heavily based on the framework described by Mokhov et al. in their prize-winning papers[@mokhov2018build;@mokhov2020build].  I outline this framework below.

## Key Concepts


### Asset

An asset is any file or directory consumed or produced by the build.  In the context of this project, examples of assets include:


- GWAS summary statistics.
- Reference RNAseq data from the [GTEx project.](../Data_Sources/GTEx_Project/GTEx_RNAseq_Data.md)
- Plots illustrating the results of applying [MAGMA](../Bioinformatics_Concepts/MAGMA_Overview.md) to GWAS summary statistics.
- Tables of category heritability weights $\{\tau_k\}$ and their associated p-values produced by [Stratified Linkage Disequilibrium Score Regression](../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md).


### Task

A task is an operation that constructs a particular asset.

Here is the source code for the `Task` base class:

```python title="Task Class Definition"
--8<-- "mecfs_bio/build_system/task/base_task.py"
```


- The `meta` property returns a metadata object (see below) describing the asset created by the Task. This metadata object must include a string `AssetId` uniquely identifying the asset.
- The `deps` property returns a list of pre-requisite Tasks that must be executed prior to the current Task.  The build system uses `deps` to construct the dependency graph of Tasks.
- The key method is `execute`. Subclasses must override this method to specify how to construct their assets.  The `fetch` parameter is a special callback passed to `execute` by the build system.  Instead of directly accessing its dependencies, `execute` should use `fetch` to access dependencies via the build system


Concrete Task subclass classes are defined [here][mecfs_bio.build_system.task].


### Rebuilder

Given a `Task` that generates an `Asset`, together with a data storage object called `Info`, the job of a Rebuilder is to decide whether the current version of the `Asset` is up-to-date.  If so, that `Asset` can be directly returned without executing the `Task`.  If not, the rebuilder uses the `Task` to materialize an up-to-date version of the asset.

Here is source code for the `Rebuilder` base class:

```python
--8<-- "mecfs_bio/build_system/rebuilder/base_rebuilder.py"
```


Currently, there is one concrete implementation of `Rebuilder`, called the [VerifyingTraceRebuilder][mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_rebuilder_core].  It uses file hashes to decide whether an `Asset` is up-to-date.


### Scheduler

Given on or more target assets requests by the user, it is the job of the scheduler to determine which tasks need to be run in what order to produce those assets.  The scheduler delegates the actual running of these tasks to the Rebuilder.

Currently, there is one concrete scheduler: the [topological scheduler][mecfs_bio.build_system.scheduler.topological_scheduler].  The topological scheduler constructs a directed acyclic graph of the dependencies of the requested assets, then traverses this graph in [topological order](https://en.wikipedia.org/wiki/Topological_sorting).


