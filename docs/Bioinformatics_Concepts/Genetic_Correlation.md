# Genetic Correlation

## Standard Definition

Like [heritability](Heritability.md), genetic correlation is typically defined in the context of an additive genetic model. Let $Y_A$ and $Y_B$ denote two phenotypes of interest.  Let $G_A,G_B$, and $E_A, E_B$ denote the respective genetic and environmental components of the phenotypes.

The model is

$$
\begin{align}
Y_A &= E_A + G_A,  \label{model1} \\ 
Y_B &= E_B + G_B.  \label{model2}
\end{align}
$$

In the context of this model, genetic correlation is defined as the Pearson correlation of the genetic components:

$$
\rho_g:= \mathrm{Corr}(G_A, G_B).
$$


## When Are Traits Genetically Correlated?

What does it mean biologically when two traits are genetically correlated? The most straightforward case is that of pleiotropy: an overlapping set of genetic variants drives the two phenotypes.  This corresponds to the causal graph $A	\leftarrow g\rightarrow B$, where $g$ denotes the shared variants.  Mediated pleiotropy is also possible, in which variants cause one trait, which causes the other: $g\rightarrow A\rightarrow B$.

Besides these straightforward cases, there are more exotic possible causes of genetic correlation, as discussed [here](https://gcbias.org/2016/04/19/what-is-genetic-correlation/).  Briefly,

- Two traits can be genetically correlated because genetics affects the behavior of a parent, which affects the phenotype of their child.
- Two traits can be genetically correlated because individuals with these traits tend to mate at a higher rate than would be expected under random mating.

## Genetic Covariance

Some applications require the calculation of the genetic covariance between two traits.  In the context of the model of $(\ref{model1},\ref{model2})$, the genetic covariance is $\mathbb{Cov}(G_A, G_B)$.  Note that genetic covariance depends strongly on how the traits are scaled.