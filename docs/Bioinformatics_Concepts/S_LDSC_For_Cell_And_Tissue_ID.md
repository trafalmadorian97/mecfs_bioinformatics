---
hide:
- navigation
- toc
---
# S-LDSC

Stratified Linkage Disequilibrium Score Regression (S-LDSC)[@finucane2018heritability;@finucane2015partitioning] is an extension of Linkage Disequilibrium Score Regression [(LDSC)](LDSC.md)[@bulik2015ld]. S-LDSC can be used to generate hypotheses about the key cells and tissues underlying a phenotype.



## High-level summary


### Review of LDSC

Recall that LDSC assumes isotropic polygenicity: heritability is evenly distributed across the genome. Mathematically this means that in LDSC  for all genetic variants $i$ and $j$:

$$
\begin{align}
\mathbb{Var}(\beta_i)&=Q\\
\mathbb{Cov}(\beta_i,\beta_j) &= 0 & \text{ when }i\ne j
\end{align}
$$

where $\beta_i$ is the true effect of variant $i$ on the phenotype and $Q\in\mathbb{R}$ is a constant.


### Overview of S-LDSC

S-LSDC instead assumes that the genome is divided into regions, and isotropic polygenicity holds within each region. Mathematically, this means that S-LDSC for all genetic variants $i$ and $j$:

$$
\begin{align}
\mathbb{Var}(\beta_i)&=\sum_{k} \tau_k I_{i \in C_k}\\
\mathbb{Cov}(\beta_i, \beta_j) &= 0 & \text{ when }i\ne j
\end{align}
$$

Where:


- The $\{C_k\}$ are sub of genetic variants.  The $C_k$ can describe functional chromosomal regions, like promoters, enhancers, or other regulatory regions.  They can also reflect gene-tissue expression.  For instance, a $C_k$ could mark SNPs that are near genes that have been observed to be over-expressed in the liver according to a [GTEx RNAseq dataset](../Data_Sources/GTEx_RNAseq_Data.md). The $C_k$ may overlap.
- The $\{\tau_k\}$ are heritability weights. $\tau_k\in\mathbb{R}$ measures the effect on the GWAS signal of a genetic variant belonging to category $C_k$.

## Procedure

### Review of LDSC Procedure

Recall that LDSC runs a univariate regression based on the approximate equation:

$$
 \mathbb{E} (\chi^2_i) \approx N\frac{h^2}{M}l_N+1,
$$

where

- $\chi^2_i$ is the $\chi^2$ statistic of the maginal GWAS regression for genetic variant $i$.
- $M$ is the number of genetic variants in the GWAS.
- $N$ is the number of individuals in the GWAS.
- $l_i:=\sum_k r_{ik}^2$ is the linkage disequilibrium score of genetic variant $i$, where $r_{ik}$ denotes the correlation between variants $i$ and $k$.

This regression provides an estimate of $h^2$, the trait's [heritability](Heritability.md).



### S-LDSC Procedure

Analogously, S-LDSC runs a multivariate regression based on the approximate equation:

$$
\mathbb{E} (\chi^2_i) \approx   N \sum_k \tau_k l_{i,k}+1,
$$
 
where the notation is the same as above, with the addition that

- $l_{i,k}:=\sum_j r_{i,j}^2 I_{j\in C_k}$ is the linkage disequilibrium score of genetic variant $i$ restricted to category $k$.

This regression provides estimates of the heritability weights $\{\tau_k\}$.

### The output of S-LDSC

- Using the $\tau_k$ weights, we can estimate the overall heritability as $h^2=\sum_i\sum_k \tau_k I_{i\in C_k}$.
- We can evaluate the proportion of heritability due to category $k'$ as $\frac{\sum_i \tau_{k'} I_{i\in C_{k'}}}{\sum_i\sum_k \tau_k I_{i\in C_k}}$.  We can compare this to the proportion of genetic variants in category $k$, which is $|C_k|/M$.  This comparison gives us an estimate of the enrichment of heritability in the category.
- We can run statistical tests to evaluate the evidence that $\tau_k\ne 0$.  If this evidence is strong, we can argue that category $k$ reflects a meaningful grouping of genetic variants that differentially affect the phenotype.
 

