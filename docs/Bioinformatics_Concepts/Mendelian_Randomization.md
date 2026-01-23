

# Mendelian Randomization

[//]: # (To discuss: Formalism of Causal inferece.  How the three instrumental variables assumptions dont let us do causal inference without a `parametric` model.  Wald ratio.  Wald ratio error.  Other methods.)
## Categories of Data Science
Hernan et al.[@hernan2019second] assert that the core activities of the data scientist are, in increasing order of difficulty:


1. **Description**: Computation of summary statistics, to enable large, complex datasets to be easily comprehended.
2. **Prediction**: Estimating the distribution of one variable conditional on another. Many applications of machine learning are prediction problems.  An example is the calculation of a person's odds of heart disease, conditional on their demographics, genetics, and [LDL](../Analysis/Verma_et_al_(LDL)/LDL_Overview.md) level.
3. **Causal Inference**:   Reasoning about counterfactuals.   For example: "If all patients with high LDL were put on statins, occurrence of heart disease would be reduced by $X$ percent".  



[//]: # (Since causal effects are the core of scientific knowledge, there is understandable interest in specifying the criteria according to which valid causal inferences can be made.)

## Methods of causal inference

The randomized controlled trial (RCT) is the gold standard causal inference method because it allows causal inference under minimal assumptions.

Unfortunately, RCTs are costly and time-consuming.  This motivates alternative causal inference techniques, which are cheap and fast, but require additional assumptions.  Mendelian Randomization (MR) is such a technique.

[//]: # (Unfortunately, there are a large number of important causal questions for which RCTs are not viable.  To address  these questions, )


## Mendelian Randomization

As an illustrative example, consider the question of how LDL affects heart disease. We might notice from epidemiological studies that people with heart disease tend to have higher LDL than people without it. We recall, however, that correlation is not causation. Thus we apply MR, in hopes of determining whether there is a true causal effect.
w

### The three main MR assumptions

The figure below from Hartley et al.[@hartley2022guide] summarizes MR:

![mr-diagram](https://github.com/user-attachments/assets/afcab2e5-b9d0-423f-a43b-0b4f47b172e4)

We apply MR to estimate the causal effect of an exposure (e.g. LDL) on an outcome (e.g. heart disease).  MR requires a genetic variant ("Genetic Instrument") associated with the exposure.

MR requires three main assumptions:

- **IV1**: The instrument must be associated with the exposure.
- **IV2**: The instrument cannot share a cause with the exposure or the outcome.
- **IV3**: The instrument cannot be causally associated with the outcome, except via the exposure.



Roughly speaking, in MR, natural variation in the instrument plays a role analogous to the random treatment assignment in an RCT.  



A key advantage of MR  is that the presence of confounders of the exposure and outcome does not affect the validity of the inference. 


### The Fourth MR Assumption

IV1, IV2, and IV3 are not sufficient to uniquely determine the causal effect of the exposure on the outcome (See Hernan and Robins Chapter 16[@hernan2010causal]).  Another assumption is required.  There are numerous possible variants forms of this forth assumption, which can largely be grouped into 

- Homogeneity assumptions: roughly, the causal effect of the exposure on the outcome is the same for all individuals.
- Monotonicity assumptions: roughly, the instrument moves the exposure in the same direction for all individuals.

