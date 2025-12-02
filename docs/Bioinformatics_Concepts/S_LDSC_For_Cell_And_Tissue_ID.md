# Stratified Linkage Disequilibrium Score Regression 

Stratified Linkage Disequilibrium Score Regression (S-LDSC)[@finucane2018heritability;@finucane2015partitioning] is an extension of [LDSC](LDSC.md).  S-LDSC relaxes some of the assumptions of LDSC, making it more broadly applicable.  Unlike LDSC, S-LDSC can also be used to generate hypothesis about the key cell-types and tissues underlying a phenotype.

## High-level summary

Recall that LDSC assumes isotropic pleiotropy: heritability is evenly distributed across the genome.

Mathematically, in LDSC we have

$$
\mathbb{Var}(\beta_i)=Q 
$$

where $\beta_i$ is the true GWAS effect size and $Q$ is some constant.




S-LSDC instead assumes that piece-wise isotropic pleiotropy: the genome is divided into a number of regions, and heritability is evenly distributed within each region.  However, heritability may be much more strongly enriched in one region than another.

Mathematically, in S-LDSC we have:

$$
\mathbb{Var}(\beta_i)=\sum_{k} \tau_k I_{i \in C_k}
$$

Where:
- The $\{C_k\}_{k}$ are subsets of the SNPs.  The $C_k$ can describe functional chromosomal regions, like promoters, enhancers, or other regulatory regions.  They can also reflect the outcome of gene-tissue expression studies.  For instance, a $C_k$ could mark SNPs that are near genes that have been observed to be over-expressed in the liver according to a [GTEx RNAseq dataset.](../Data_Sources/GTEx_project).
- $\tau_k$ is a weight reflecting the enrichment of GWAS signal due to the category $C_k$.

## Procedure

Recall that LDSC run a univariate regression using the equation

$$
 \mathbb{E} (\chi^2_j) \approx \frac{h^2}{M}l_jN+1
$$

to estimate $h^2$, the trait's SNP-heritability.

Analogously, S-LDSC runs a multivariate regression using the equation

$$
\mathbb{E} (\chi^2_j) \approx 1+ N + N \sum_k \tau_k l_{j,k}
$$
 

