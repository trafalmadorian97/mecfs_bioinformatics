# Cross Trait Linkage Disequilibrium Score Regression
Cross Trait Linkage Disequilibrium Score Regression (CT-LDSC)[@bulik2015atlas] is an [extension of Linkage Disequilibrium Score Regression](LDSC.md) (LDSC) that allows the estimation of [genetic correlation](Genetic_Correlation.md).

## Derivation

### Data Generating Model

Recall that LDSC assumes a linear data-generating model for the trait of interest.  Similarly, CT-LDSC assumes two linear data-generating models.

$$
\begin{align}
y_1&= Y\beta + \delta
y_2&= Z\gamma + \epsilon
\end{align}
$$


where:

- There are $M  \gg 0$ genetic variants.  
- There are $N_1 \gg 0$ individuals in the GWAS of trait 1 and $N_2 \gg 0$ individuals in the GWAS of trait 2.
- $y_1\in \mathbb{R}^{N_1}$ and $y_2 \in \mathbb{R}^{N_2}$ are the vectors of phenotypes for trait 1 and trait 2 respectively.
- $Y\in\mathbb{R}^{N_1\times M}$ and $Z\in\mathbb{R}^{N_2\times M}$ are the genotype matrices from the two GWAS, normalizd to have columns with mean 0 and variance 1.
- $\beta,\gamma\in\mathbb{R}^M$ are the vectors of true genetic effect sizes for each genetic variant for the two traits. 
- $Y\beta\in\mathbb{R}^{N_1}$ and $Z\gamma\in\mathbb{R}^{N_2}$ are thus the vectors of true genetic effects for all individual in the two GWAS.
- $\epsilon\in\mathbb{R}^{N_1} and $\delta\in\mathbb{R}^{N_2}$ are the vectors of non-genetic effects for all individuals in the two GWAS.
