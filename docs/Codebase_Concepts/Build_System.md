# Build System

## Motivation
In data science, the final product of an analysis is often the result of long pipeline of many steps.  These steps can be a mixture of pure data cleaning operations with more complex analysis and modeling.  Moreover, when the data size is nontrivial, the steps and therefore the whole pipeline can be quite slow.  This scenario creates two challenges:

- **Iteration**: It is rare to run a pipeline once, produce an analysis, and be done with it. Usually, it is necessary to repeatedly tweak the steps, re-run the pipeline, and examine the result.  To avoid wasting time, it is therefore desirable that after each change, only the steps impacted by the change should be rerun.
- **Lineage**:  Given the complexity of many data science workflows, there is considerable room of error.  It is therefore desirable to be able to interrogate the final product of a workflow to trace its "lineage": the precise series of steps that produced it.

These challenges motivate the development of a data science build system.

## Framework 

The build system used in this project is heavily based on the framework described by Mokhov et al. in their prize-winning papers[@mokhov2018build] [@mokhov2020build].  

## Key Concepts

### Asset

An asset is any file or directory relevant to the build.  Examples of assets include:


- GWAS summary statistics
- Reference RNAseq data from the GTEx project
- Plots illustrating the results of applying MAGMA to GWAS summary statistics.


### Task

A task is an operation that constructs a particular asset.

Here is the source code for the `Task` base class:

```python title="Task Class definition"
--8<-- "mecfs_bio/build_system/task/base_task.py"
```


- The `meta` property returns a metadata object (see below) describing the asset created by the Task. This metadata object must include a string `AssetId` uniquely identifying the asset.
- The `deps` property returns a list of pre-requisite Tasks that must be executed prior to the current Task.  The build system uses `deps` to construct the dependency graph of Tasks.
- The key method is `execute`. Subclasses must override this method to specify how to construct their assets.  The `fetch` parameter is a special callback passed to `execute` by the build system.  Instead of directly accessing its dependencies, `execute` should use `fetch` to access dependencies via the build system


Concrete Task subclass classes are defined [here][mecfs_bio.build_system.task].


