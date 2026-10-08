# Liability Threshold Model

## Overview

Modern statistical genetics treats the genetic contribution to common disease as quantitative and continuous. In contrast, many common disease phenotypes are binary: a person either has the disease or does not.  How can we bridge this gap? The liability threshold model is simple, elegant solution.


Under the liability threshold model, the genetic contribution to disease $G$ combines additively with the environmental contribution $E$ to produce the latent disease liability $L$. When $L$ exceeds the threshold $\tau$, the patient expresses the disease phenotype $Y=1$.  Otherwise, the patient expresses the normal phenotype $Y=0$.  In equations, we have


$$
\begin{align}
L&=G+E\\
Y&= 1_{L>\tau}.
\end{align}
$$


It is typical to model $G$ and $E$ as independent normal random variables.  This allows the liability threshold model to analyzed via the theory of the [truncated normal distribution](https://en.wikipedia.org/wiki/Truncated_normal_distribution).

The diagram below illustrates a typical liability threshold model.


{{
static_img_embed("docs/_figs/liability_threshold_model_diagram.svg", width="100%")
}}



A key advantage of the liability threshold model is that it allows many statistical techniques originally developed for quantitative phenotypes to be applied to the binary phenotypes. One just applies such a technique to the underlying liability $L$ instead of the observed phenotype $Y$.


## Liability-scale heritability

Where we model a binary trait as being governed by the liability threshold model, it is frequently of interest to compute [heritability](Heritability.md) on the liability scale.  This is defined as:


$$
\begin{align}
h^2_L:=\frac{\mathrm{Var}(G)}{\mathrm{Var}(L)}.
\end{align}
$$


Thus, we simply apply the standard definition of heritability to the liability $L$ instead of the phenotype $Y$.

## Ascertained case-control studies

For practical reasons, often participants in binary-trait GWAS are selected not by sampling the general population uniformly at random, but through special ascertainment procedures that differ for cases and for controls.   Cases are recruited from doctor's offices, hospitals, patient advocacy groups, etc., while controls are selected from large scale biobanks so that they are good match for controls in their ancestry, age, gender, etc[^matching].  

The above-described ascertainment interacts with the liability threshold model.  If we assume that liability $L$ is normally distributed in the general population, it will not be normally distributed in the ascertained GWAS sample



[^matching]: This matching reduces, but does not eliminate, the need to control for confounders.