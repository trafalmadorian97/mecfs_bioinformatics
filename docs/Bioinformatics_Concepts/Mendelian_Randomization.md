

# Mendelian Randomization

[//]: # (To discuss: Formalism of Causal inferece.  How the three instrumental variables assumptions dont let us do causal inference without a `parametric` model.  Wald ratio.  Wald ratio error.  Other methods.)
## Categories of Data Science
Hernan et al.[@hernan2019second] assert that the core activities of the data scientist are, in increasing order of difficulty:


1. **Description**: Computation of summary statistics, to enable large, complex datasets to be easily comprehended.
2. **Prediction**: Estimating the distribution of one variable conditional on another. Many of the important problems to which machine learning has been applied are prediction problems.  An example would be the classification of patients as having or not having diabetic retinopathy based on images of their retinas.
3. **Causal Inference**:   Reasoning about counterfactuals.   For example: "If all newly-diagnosed diabetes patients were prescribed metformin, occurrence of diabetic retinopathy would be reduced by $X$ percent".  Causal inference allows us to understand when we can make statements about the consequences of an intervention. 



[//]: # (Since causal effects are the core of scientific knowledge, there is understandable interest in specifying the criteria according to which valid causal inferences can be made.)

## Methods of causal inference

The gold standard method for causal inference is the randomized controlled trial (RCT).  RCTs are rightly privileged for allowing causal conclusions to be drawn while requiring minimal assumptions.

Unfortunately, there are a large number of causal questions for which RCTs are not viable.  For these questions, it is necessary to apply other causal inference techniques, which require additional assumptions.  Mendelian Randomization (MR) is such a technique.

## Mendelian Randomization

As an illustrative example of MR, consider the question of determining the causal effect of [Low Density Lipoprotein (LDL)](../Analysis/Verma_et_al_(LDL)/LDL_Overview.md) on heart disease. We might notice from epidemiological studies that people who experience heart disease tend to have higher LDL than those who do not. We recall, however, that correlation is not causation. Thus we apply MR, in hopes of determining whether there is a true causal effect.

### The three main MR assumptions

The figure below from Hartley et al.[@hartley2022guide] summarizes Mendelian Randomization (MR):

![mr-diagram](https://github.com/user-attachments/assets/afcab2e5-b9d0-423f-a43b-0b4f47b172e4)

We apply MR to estimate the causal effect of an exposure (e.g. LDL) on an outcome (e.g. heart disease).  MR requires a genetic variant (or "Genetic Instrument") that must be associated with the exposure (**IV1** in the diagram).

MR requires additional assumptions:

- **IV2**: The instrument can share a cause with the exposure or the outcome.
- **IV3**: The instrument cannot be causally associated with the outcome, except via the exposure.

A key advantage of MR compared to other causal inference techniques is that the presence of confounders of the exposure and outcome does not affect the validity of the inference. 


Roughly speaking, natural variations in the genetic instrument affect the exposure, which affect the outcome.  By observing individuals with these natural variations, we estimate the effect of the exposure on the outcome.

Thus, variation in the instrument plays a role analogous to random treatment assignment in an RCT.  


### The Fourth MR Assumption




