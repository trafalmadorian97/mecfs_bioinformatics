

# Mendelian Randomization

[//]: # (To discuss: Formalism of Causal inferece.  How the three instrumental variables assumptions dont let us do causal inference without a `parametric` model.  Wald ratio.  Wald ratio error.  Other methods.)
## Categories of Data Science
Hernan et al.[@hernan2019second] assert that the core activities of the data scientist are, in increasing order of difficulty:
1. **Description**: The computation of summary statistics, which enables key aspects of large, complex datasets to be easily comprehended. Given the abundance of data in today's world, this activity is fundamental.
2. **Prediction**: The fitting of models that estimate the distribution of one variable conditional on another. Many of the important problems to which machine learning has been applied require only prediction.  An example would be the classification of patients as having or not having diabetic retinopathy based on images of their retinas. It has also been argued that certain kinds of important policy or economic decisions can be made using only prediction[@kleinberg2018human].
3. **Causal Inference**:   Reasoning about counterfactuals. Causal inference allows us to understand when we can make statements about the consequences of an intervention.   For example: "If all newly-diagnosed diabetes patients were prescribed metformin, occurrence of diabetic retinopathy would be reduced by $X$ percent".  If we can provide strong evidence supporting such a statement, we are said to have demonstrated a causal effect.

Since causal effects are the core of scientific knowledge, there is understandable interest in specifying the criteria according to which causal inference is valid.

## Methods of causal inference

The gold standard for causal inference is the randomized controlled trial (RCT).  RCTs are rightly privileged for allowing causal conclusions to be drawn from minimal assumptions.

Unfortunately, there are a large number of causal questions for which RCTs are not viable.  For these questions, it is necessary to apply other causal inference techniques, which depend on additional assumptions.  Instrumental variables is an example of such a technique.

## Instrumental variables




