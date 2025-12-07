# Build System

## Motivation
In data science, the final product of an analysis is often the result of long pipeline of many steps.  These steps can be a mixture of pure data cleaning operations with more complex analysis and modeling.  Moreover, when the data size is nontrivial the steps and therefore the whole pipeline can be quite slow.  This scenario creates two challenges:

- **Iteration**: It is rare to run a pipeline once, produce an analysis, and be done with it. Usually, it is necessary to repeatedly tweak the steps, re-run the pipeline, and examine the result.  To avoid wasting time, it is therefore desirable that after each change, only the steps impacted by the change should be rerun.
- **Lineage**:  Given the complexity of many data science workflows, there is considerable room of error.  It is therefore desirable to be able to interrogate the final product of a workflow to trace its "lineage": the precise series of steps that produced it.

These challenges motivate the development of a build system.