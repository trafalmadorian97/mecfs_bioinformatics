# Genetic Correlation

## Standard Definition

Like [heritability](Heritability.md), genetic correlation is typically defined in the context of an additive genetic model. Let $Y_A$ and $Y_B$ denote the two phenotypes of interest.  Let $G_A,G_B$, and $E_A, E_B$ denote the respective genetic and environmental components of the phenotypes.

The model is

$$
\begin{align}
Y_A &= E_A + G_A,\\
Y_B &= E_B + G_B.
\end{align}
$$

In the context of this model, genetic correlation is defined as the Pearson correlation of the genetic components of the two phenotypes:

$$
\rho_g:= \mathrm{Corr}(G_A, G_B)
$$


## When Are Traits Genetically Correlated?

What does it mean biologically when two traits are genetically correlated? The most straightforward case is that of pleiotropy: the overlapping set of genetic variants are driving the two phenotypes.  This corresponds to the causal graph A<-g->B where g denotes the shared variants.  Mediated pleiotropy is also possible, in which genetics causes one trait which causes another: g->A->B.

Beside these straightforward cases, there are also a number of more exotic possibile causes of genetic correlation, as discussed [here](https://gcbias.org/2016/04/19/what-is-genetic-correlation/).  Briefly

- Two traits can be genetically correlated because genetics affects the behavior of a parent, which affects the phenotype of their child.
- Two traits can be genetically correlated because individuals with these traits tend to have offspring at a higher rate than would be expected under random mating.