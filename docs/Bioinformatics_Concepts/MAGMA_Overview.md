# MAGMA
## MAGMA Overview
[MAGMA](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219) is a tool for the analysis of GWAS results. MAGMA can operate on GWAS summary statistics and does not require raw GWAS data.


[MAGMA](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219) performs two main steps:

1. Gene analysis,
2. Gene property analysis.


## Gene Analysis
### Purpose
GWAS summary statistics consist of [variant](https://en.wikipedia.org/wiki/Genetic_variant)-level regression coefficients and standard errors. The purpose of the MAGMA gene-analysis step is to convert these variant-level statistics to gene-level statistics, so we can judge which genes are likely to affect the phenotype.

### Requirements
Gene analysis requires:

1.  GWAS summary statistics,
2. A rule associating variants with genes (often just proximity),
3. Data describing the [linkage disequilibrium](https://en.wikipedia.org/wiki/Linkage_disequilibrium) structure of the variants studied by the GWAS.

### Mathematical Overview


 Let $Y$ be the phenotype of interest in a GWAS.  Let $X_i$ be the ith variant.


In a GWAS, we estimate

$$
Y= \beta_i X_i +\epsilon_i,
$$

for all $i$, where $\beta_i$ is the regression coefficient of $Y$ on $X_i$, and $\epsilon_i$ is the regression error. 


Let $G$ be the set of variants that regulate or modify the gene of interest.  We wish to measure the strength of the evidence that $\beta_i\ne 0$ for some $i\in G$. That is, we wish to judge how likely it is that the gene affects the phenotype.


To evaluate the evidence, [MAGMA uses as its test statistic](https://vu.data.surfsara.nl/s/VeuWKUwd0rz6AZD?dir=/&editing=false&openfile=true)[@MagmaSNPWiseNote;@de2015magma] $\sum_{i\in G}Z_i^2$, where $Z_i$ is the Z-statistic of variant $i$'s GWAS regression. 

This test statistic has a [generalized chi-squared distribution](https://en.wikipedia.org/wiki/Generalized_chi-squared_distribution) under the null hypothesis.  The details of its distribution depend on the correlations between the individual $Z$-statistics, which in turn depends on the linkage disequilibrium structure of the variants under study.  This is why MAGMA gene analysis requires linkage disequilibrium reference data.

MAGMA converts the test-statistic to a p value via a[ numerical-integration procedure](https://vu.data.surfsara.nl/s/VeuWKUwd0rz6AZD?dir=/&editing=false&openfile=true). A small p value indicates strong evidence that the gene affects the phenotype.

## Gene Property Analysis

### Purpose

While knowing which genes are associated with the phenotype is of independent interest, understanding a phenotype requires knowing the biological systems that underlie it. To facilitate this, MAGMA provides a gene property module, which aims to answer the question: "are genes 
associated with the phenotype over-represented in particular biological systems?"

## Requirements
The gene property analysis module requires:

1. The output of the gene analysis step, assigning p values to genes.
2. A dataset of gene properties.  This dataset measures the extent to which different genes participate in various biological systems.  For example, tissue-specific RNAseq data from [GTEx](../Data_Sources/GTEx_RNAseq_Data.md) could be used.  A high quantity of RNA transcripts of a gene in a given tissue shows that the gene plays a role in that tissue.
 
## Mathematical Overview
- Let $Z_i$ be the $z$-score of the association of gene $i$ with the phenotype.
- Let $E_{i,j}$ be the measure of the participation of gene $i$ in biological system $j$.
- Let $m$ be a vector of control covariates. 

For each biological system $j$, MAGMA fits the regression

$$
\begin{align}
Z_i = \beta_0 + \beta_{j} E_{i,j} + (\beta^m_{j})^T m + \epsilon_{i,j}, \label{magma_gene_prop}
\end{align}
$$

where the $\beta$ are regression coefficients and $\epsilon_{i,j}$ is the regression error.

The null hypothesis $\beta_{j}= 0$ is tested.

Rejecting this null would indicate that the degree to which the gene $i$ participates in the biological system $j$ is predictive of the extent of association of the gene with the phenotype, $Z_i$.  This would suggest that the biological system is related to the phenotype.

## Limitations

While MAGMA is a very useful approach for detecting the key biological systems involved in a disease or trait, it has a number of limitations.  These include:

1.  It is fundamentally limited by the availability and accuracy of data associating genes with biological systems.  For instance, if tissue-specific RNA-seq data is used, but RNAseq has not been performed on the key tissue involved in the disease, MAGMA will not produce useful results.
2. It is unknown how to translate raw reference data, derived from RNAseq or another method, into $E_{i,j}$ values to optimize the resulting MAGMA signal.  This question is its own research topic[@li2025benchmarking].
3. Relatedly, MAGMA relies on the assumption that if a biological system is important to a phenotype, the degree of participation of a gene in that biological system should be linearly related to the association of the gene with the phenotype.  While this is intuitively reasonable, it does not necessarily hold: biology is inherently nonlinear.
4. Since MAGMA in its original form mostly uses proximity to associate genes with variants, it cannot incorporate long-distance regulatory effects.  If such effects are important to a disease, MAGMA could give misleading results.


## References
The above discussion is based on:

[//]: # ( [Generalized chi-squared distribution]&#40;https://en.wikipedia.org/wiki/Generalized_chi-squared_distribution&#41;)

[FUMA Tutorial](https://fuma.ctglab.nl/tutorial#snp2gene)

[MAGMA Paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219) 

[MAGMA Site](https://cncr.nl/research/magma/)

[Note on Magma SNP-wise model](https://vu.data.surfsara.nl/s/VeuWKUwd0rz6AZD?dir=/&editing=false&openfile=true)

[Paper Critiquing the original MAGMA SNP-wise model](https://web.archive.org/web/20240724122326/https://www.biorxiv.org/content/10.1101/2020.08.20.260224v1) (MAGMA was updated in response to this critique.)


\bibliography
