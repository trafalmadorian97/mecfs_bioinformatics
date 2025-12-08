# Stratified Linkage Disequilibrium Score Regression 

Stratified Linkage Disequilibrium Score Regression (S-LDSC)[@finucane2018heritability;@finucane2015partitioning] is an extension of Linkage Disequilibrium Score Regression [(LDSC)](LDSC.md).  S-LDSC relaxes some of the assumptions of LDSC, making it more broadly applicable. S-LDSC can be used to generate hypotheses  about the key cell-types and tissues underlying a phenotype.

## High-level summary

### Review of LDSC

Recall that LDSC assumes isotropic pleiotropy: heritability is evenly distributed across the genome.

Mathematically, in LDSC we have that for all genetic variants $i$

$$
\mathbb{Var}(\beta_i)=Q 
$$

where $\beta_i$ is the true effect size and $Q$ is some constant.


### Overview of S-LDSC

S-LSDC instead assumes that genome is divided into regions, and isotropic pleiotropy holds  within each region.

Mathematically, in S-LDSC we have that for all genetic variants $i$:

$$
\mathbb{Var}(\beta_i)=\sum_{k} \tau_k I_{i \in C_k}
$$

Where:


- The $\{C_k\}_{k}$ are subsets of the SNPs.  The $C_k$ can describe functional chromosomal regions, like promoters, enhancers, or other regulatory regions.  They can also reflect gene-tissue expression.  For instance, a $C_k$ could mark SNPs that are near genes that have been observed to be over-expressed in the liver according to a [GTEx RNAseq dataset.](../Data_Sources/GTEx_Project/GTEx_RNAseq_Data.md).
- The $\{\tau_k\}$ are weights, where $\tau_k$ measures the effect on the GWAS signal of category $C_k$.

## Procedure

### Review of LDSC Procedure

Recall that LDSC runs a univariate regression using the equation

$$
 \mathbb{E} (\chi^2_i) \approx \frac{h^2}{M}l_iN+1
$$

where

- $\chi^2_i$ is the $\chi^2$ statistic of the GWAS regression of the $i$th genetic variant.
- $M$ is the number of genetic variants under consideration
- $N$ is the number of individuals in the GWAS.
- $l_i=\sum_k r_{ik}^2$ is the linkage disequilibrium score of genetic variant $i$, where $r_{ik}$ denotes the correlation between variants $i$ and $k$

This regression provides as estimate of $h^2$, the trait's [heritability](Heritability.md).



### S-LDSC Procedure

Analogously, S-LDSC runs a multivariate regression using the equation

$$
\mathbb{E} (\chi^2_i) \approx 1+ N + N \sum_k \tau_k l_{i,k}
$$
 
where the notation is the same as above, and where

- $l_{i,k}:=\sum_j r_{i,j}^2 I_{j\in C_k}$ is the linkage disequilibrium score of genetic variant $i$ with category $k$.

The regression provides estimates of the category heritabilities weights $\{\tau_k\}$.

Using these weights, one can estimate the overall heritability as $h^2=\sum_i\sum_k \tau_k I_{i\in C_k}$.

