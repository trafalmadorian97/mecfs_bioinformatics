---
hide:
- navigation
- toc
---
# Linkage Disequilibrium Score Regression 
Linkage Disequilibrium Score Regression[@bulik2015ld] (LDSC) is a technique for estimating  [heritability](Heritability.md) from GWAS summary statistics.  LDSC is ubiquitous, but its usefulness depends on the validity certain modeling assumptions. This page includes both a high-level summary and a detailed derivation of LDSC.

## Intuition and High-level summary

The core idea of LDSC is illustrated below:

![ldsc-schematic](https://github.com/user-attachments/assets/7a642f44-2fb9-4647-888a-23c1a344d12d)

For each SNP, we compute a quantity called Linkage-Disequilibrium Score, or LD score, which measures the strength of the correlation between this SNP and other SNPs.  Two key insights are central to LDSC:

- SNPs with higher LD scores will tend to be more strongly associated with the GWAS phenotype, and thus will have higher [Wald test statistics](https://en.wikipedia.org/wiki/Wald_test).  This follows because if a SNP is not itself causal, it can be associated with the GWAS phenotype through correlation with a causal SNP.  A SNP with a higher LD score has a higher chance of being correlated with a causal SNP than a SNP with a lower LD score.
- The strength of the relationship between LD score and the Wald test statistic depends on heritability.  A more heritable trait (right panel) will exhibit a stronger relationship than a less heritable trait (left panel).  This follows because for a more heritable trait, SNPs have a greater influence over the phenotype, so that increasing the extent to which a SNP is in LD with other SNPs will on average increase its association with the phenotype by a larger amount.


Using these two insights, we can estimate heritability from the slope of regression line of Wald test statistics on LD scores.

[//]: # (### Core Assumption)

[//]: # (Roughly speaking, the core assumption of LDSC is that:)

[//]: # ()
[//]: # (*Heritability is evenly distributed across the genome*)

[//]: # ()
[//]: # (The accuracy of LDSC will depend on how close this assumption is to being satisfied)

[//]: # ()
[//]: # (### Core Insight)

[//]: # ()
[//]: # (The core insight of LDSC is that)

[//]: # ()
[//]: # (*A SNP in strong linkage disequilibrium with other SNPs is more likely to be significantly associated with the phenotype than a SNP in weak linkage disequilibrium*)

[//]: # ()
[//]: # (This follows because in a GWAS, a SNP can be associated with the phenotype either)

[//]: # ()
[//]: # (-  by being causal, or)

[//]: # (-  by being in linkage disequilibrium with a causal SNP.)

[//]: # ()
[//]: # (Because we have assumed heritability is evenly distributed, a SNP in strong LD has a higher chance of being in LD with the causal SNP than a SNP in weak LD.)

[//]: # ()
[//]: # (LDSC uses this core insight to estimate heritability.)


[//]: # (### How it works)

[//]: # ()
[//]: # (The relationship between strength of LD and significance of association is determined by the heritability of the trait.  To see this intuitively, note that if the trait is not heritable, all SNPs will have negligible association with the phenotype.  In contrast, if the trait is highly heritable, the gap in strength of association between SNPs in strong LD and SNPs in weak LD will be large, since SNPs in strong LD will be likely to be in LD with many causal SNPs.  LDSC uses this principle of the relationship between LD and significance to go backward from the pattern of SNP significance in a GWAS to the heritability of the trait.)

[//]: # ()
[//]: # ()







## Detailed Derivation of Method
### Data-generating model

LDSC assumes the following data generating equation:

$$
\phi = X\beta + \epsilon
$$

where:


- There are $M  \gg 0$ genotypes.  
- There are $N  \gg 0$ individuals.  
- $\phi\in\mathbb{R}^N$ is the vector of phenotypes.
- $X\in\mathbb{R}^{N\times M}$ is the genotype matrix, normalized to have columns with sample mean 0 and sample variance 1.
- $\beta\in\mathbb{R}^M$ is the vector of true SNP effect sizes.
- $X\beta \in\mathbb{R}^N$ is thus the vector of total genetic effects for each individual.
- $\epsilon\in\mathbb{R}^N$ is the vector non-genetic effects for each individual.

Furthermore, we model $X,\beta,\epsilon$ as random variables with the following properties:

$$
\begin{align}
\mathbb{Var}(\epsilon)&= (1-h^2)I  & \text{For some $h>0$} \\
\mathbb{Var}(\beta)&=\frac{h^2}{M}I \label{betavar} \\
\mathbb{E}(\epsilon)&=0\\
\mathbb{E}(\beta)&=0 \\
\mathbb{E} X &=0 \\ 
\end{align}
$$

- The rows of $X$ are independent and identically distributed.
- The distributions of columns of $X$ are "not too similar".
- The SNP $X_{i,j}$ may be highly correlated with a few other SNPs, but is uncorrelated with most SNPs.
- $\beta,\epsilon,X$ are all mutually independent.


Furthermore, define the following quantities related to Linkage Disequilibrium (LD).

- The LD between SNP $j$ and SNP $k$ is denoted by $r_{jk}:=\mathbb{E}X_{i,j}X_{i,k}$, (which does not depend on the individual $i$ by our assumption that the rows of $X$ are iid).
- The empirical LD between $j$ and $k$ is defined to be $\tilde{r}_{jk}:=\frac{1}{N}X_{:,j}^T X_{:,k}$.
- The LD score for a SNP $j$ is defined to be $l_j:= \sum_k r_{jk}^2$.



### Properties of $\hat{\beta}_j$ and $\chi^2_j$
We begin by computing the variance matrix of the genetic effects $X\beta$:


$$
\begin{align}
&\mathbb{Var}(X\beta)\\
&=\mathbb{E} X \beta \beta^T X^T\\
&= \mathbb{E} (X \mathbb{E} ( \beta \beta^T |X ) X^T)  & \text{Tower law of expectation} \\
&=\frac{h^2}{M} \mathbb{E} X X^T & \text{by (\ref{betavar})} \\
&=h^2I & \text{By independence of rows, zero mean of $X_{i,j}$}
\end{align}
$$

We can also compute the variance matrix of the phenotypes:

$$
\begin{align}
&\mathbb{Var}(\phi)\\
&= \mathbb{Var}(X\beta) + \mathbb{Var}(\epsilon) & \text{ By independence}\\
&=h^2I + (1-h^2)I\\
&= I.
\end{align}
$$

Note that the [heritability](Heritability.md) of the phenotype for an arbitrary individual $i$ is:

$$
\begin{align}
\frac{\mathbb{Var}(X\beta)_{i,i} }{\mathbb{Var}(\phi)_{i,i}}&= h^2
\end{align}
$$

justifying the use of the symbol $h^2$.

The expectation of the phenotype vector is:

$$
\begin{align}
&\mathbb{E} (\phi) \\
&= \mathbb{E}(X \beta) + \mathbb{E} (\epsilon)\\
&= 0 & \text{Since $\mathbb{E}\beta =\mathbb{E} \epsilon=0$ }
\end{align}
$$


In a GWAS, it is typical to run a separate single-variant regression for each variant. Denote by $\hat{\beta}_j$ the single-variant regression coefficient for SNP $j$.  This is given by:

$$
\begin{align}
\hat{\beta}_j &\\
&= \frac{\frac{1}{N}(\phi-\overline{\phi})^T(X_{:,j} -\overline{X_{:,j}} )}{\frac{1}{N}  (X_{:,j} -\overline{X_{:,j}} )^T (X_{:,j} -\overline{X_{:,j}} ) }    & \text{ By OLS formula}\\
&= \frac{\frac{1}{N}(\phi-\overline{\phi})^TX_{:,j}  }{  \frac{1}{N}X_{:,j}^TX_{:,j}  } & \text{Since columns of $X$ are normalized}\\
&\approx \frac{\frac{1}{N}\phi^TX_{:,j}  }{\frac{1}{N}X_{:,j}^TX_{:,j}} & \text{Since $\mathbb{E}\phi=0$ and $N$ is large}\\
&=\frac{1}{N}\phi^TX_{:,i}  & \text{Since columns of $X$ are normalized}
\end{align}
$$

The squared regression errors of the $j$th GWAS regression are given by

$$
\begin{align}
 \lVert\phi- \hat{\beta}_jX_{:,j} \rVert^2 = \lVert\phi- \frac{1}{N}(\phi^T X_{:,j}) X_{:,j} \rVert^2.
\end{align}
$$ 


We have assumed that $M$ is large, that columns of $X$ are not too similar, and that the components of $\beta$ are iid. Thus, GWAS regression on SNP $j$ will not explain any significant proportion of the variance of the phenotype.  So  


$$
\begin{align}
&\lVert\phi- \frac{1}{N}\phi^T X_{:,j} X_{:,j} \rVert^2\\
&\approx \lVert\phi\rVert^2    \\
&\approx \mathrm{trace}\left( \mathbb{Var} (\phi) \right)\\
&= N  \label{residuals}
\end{align}
$$


The standard error of the $j$th GWAS regression coefficient is





$$
\begin{align}
&\mathrm{SE}(\hat{\beta}_j) \\
&\approx \sqrt{\frac{\frac{1}{N}\lVert\phi- \frac{1}{N}\phi^T X_{:,j} X_{:,j} \rVert^2}{(X_{:,j} -\overline{X_{:,j}} )^T (X_{:,j} -\overline{X_{:,ij} )}}} & \text{Formula for OLS SE}\\
&=\sqrt{\frac{1}{N}} & (\text{\ref{residuals}) + Normalization of $X$} \label{sebeta}
\end{align}
$$






By the definition of the [Wald Test Statistic](https://en.wikipedia.org/wiki/Wald_test) for the $j$th SNP, $\chi^2_j$, we have

$$
\begin{align}
&\chi_j^2\\
&= \frac{\hat{\beta_j}^2}{\mathrm{SE}(\hat{\beta_j})^2}\\
&\approx N \hat{\beta_j} & \text{by (\ref{sebeta})} \label{wald}
\end{align}
$$

### Expectation of empirical LD scores
The expectation of the square of the empirical LD between SNPs $j$ and $k$ is:

$$
\begin{align}
&\mathbb{E} \tilde{r}_{jk}^2\\
&=\mathbb{E}\frac{1}{N^2}\sum_i \sum_q X_{i,j}X_{i,k}X_{q,j}X_{q,k}\\
&=\frac{1}{N^2}\mathbb{E}(\sum_{i\ne q} X_{i,j}X_{i,k}X_{q,j}X_{q,k} + \sum_i X_{i,j}^2X_{i,k}^2)\\
&=\frac{1}{N^2}( \sum_{i\ne q} \mathbb{E}(X_{i,j} X_{i,k })\mathbb{E} (X_{q,j X_{q,k}}) + \sum_i \mathbb{E} X_{i,j}^2 X_{i,k}^2) &\text{Independence of rows of $X$}\\
&= \frac{N-1}{N}r_{jk}^2 + \frac{1}{N}\mathbb{E}X_{1j}^2X_{1k}^2
\end{align}
$$


We write, $\mathbb{E}X_{1j}^2X_{1k}^2=1+2r_{jk}^2+\nu$ where we have used[ Isserlis's Theorem](https://en.wikipedia.org/wiki/Isserlis%27s_theorem) to compute the expectation of the product of the squares of two normal random variables, and then added error term $\nu$ to account for the non-normality of $X$.
 
Thus we have 

$$
\begin{align}
&\mathbb{E} \tilde{r}_{jk}^2 \\
&=r_{jk}^2 + \frac{1}{N}(1+r_{jk^2}^2+\nu)\\
&=r_{jk}^2+ \frac{1}{N}(2r_{jk}^2+\nu)+ \frac{1}{N}(1-r_{jk}^2)
\end{align}
$$

The expectation of the sum of the squares of the empirical LD scores with SNP j is given by

$$
\begin{align}
&\mathbb{E} \sum_k \tilde{r}_{jk}^2\\
&=\sum_{k}r^2_{jk} + \sum_{k}\frac{1}{N}(2r_{jk}^2+\nu) + \sum_k  \frac{1}{N}(1-r_{jk}^2)\\
&=l_j + \sum_{k}\frac{1}{N}(2r_{jk}^2+\nu) + \sum_k  \frac{1}{N}(1-r_{jk}^2) & \text{Definition of $l_j$}\\
\end{align}
$$

At this stage we introduce another approximation.  We have assumed that SNP $j$ has high correlation with only a small number of other SNPs.  Thus the first $O(1/N)$ term can will be much smaller than the second.  We therefor neglect this term, yielding:

$$
\mathbb{E}\sum_k \tilde{r}_{jk}^2 \approx l_j + \frac{1}{N}(M-l_j)
$$

### The LDSC Equation

With the preliminaries out of the way, now compute the expectation of the chi squared statistic of the $j$th SNP

$$
\begin{align}
&\mathbb{E} (\chi^2_j)\\
&=N \mathbb{Var}(\hat{\beta_j)} & \text{by (\ref{wald})}\\
&=N \mathbb{E} \mathbb{Var}(\hat{\beta}_j | X) + \underbrace{\mathbb{Var}\mathbb{E}( \hat{\beta}_j | X  )}_{=0}   & \text{Law of Total Variance}\\
&= N \mathbb{E}(  \mathbb{Var}(  \frac{1}{N} \phi^T X_{:,j} |X)    )  \label{varphi}   \\
&= N \mathbb{E}(  \mathbb{Var}(  \frac{1}{N} (X\beta + \epsilon)^T X_{:,j} |X)    )    \\
&= N \mathbb{E}(  \mathbb{Var}(  \frac{1}{N}( (X\beta)^TX_{:,j} + \epsilon^T X_{:,j}) |X)    )    \\
&= N \mathbb{E}(  \mathbb{Var}(  \frac{1}{N}( X_{:,j}^T (X\beta)+ X_{:,j}^T\epsilon   ) |X)    )    \\
&= N \mathbb{E}(   \frac{1}{N}( X_{:,j}^T X \mathbb{E}(\beta \beta^T|X)  X^T X_{:,j}   + X_{:,j}^T \mathbb{E}(\epsilon\epsilon^T)X_{:,j}   )     )    \\
&= N \mathbb{E} ( \frac{h^2}{MN^2} X_{:,j}^T X X^T X_{:,j}+ \frac{1}{N}  (1-h^2) )\\
&= N \mathbb{E} ( \frac{h^2}{M}\sum_k \tilde{r}^2_{jk} + \frac{1}{N}  (1-h^2) ) &\text{def of }\tilde{r}^2_{jk} \label{hstep}\\
&= \frac{h^2}{M}(Nl_j + M-l_j) + 1 - h^2\\
&= \frac{h^2}{M}l_j(N-1)+1\\
&\approx \frac{h^2}{M}l_jN+1\\
\end{align}
$$

This is the main Linkage Disequilibrium Score Regression equation.

## Intuition about derivation
The key steps in the LDSC derivation above are between equations ($\ref{varphi}$) and ($\ref{hstep}$).  These steps relate the GWAS regression coefficients to a measure of linkage disequilibrium between SNPs.  These steps are only possible because of our distributional assumptions on $\beta$ and $\epsilon$.



## Critique of assumptions

[//]: # (Numerous authors have criticised the plausibility of the assumptions underlying LDSC.  For instance ...)


We saw above that the most important of LDSC's assumptions is that  $\mathbb{Var}\beta=h^2 I$.  This can be understood as an **isotropic pleiotropy** assumption, and it amounts to the assumption that the heritability of a trait is distributed evenly across the genome, without correlation between SNPs.

How plausible is this assumption?

- On the one hand, the discovery that many traits are highly pleiotropic has been one of the most important findings of the GWAS era.  So, in a rough sense, assuming that the heritability of a trait is distributed across the genome is not unreasonable.
- On the other hand, the assumption of perfectly uniform pleiotropy strains plausibility.  For most traits, heritability is concentrated in certain key regions.  In autoimmune diseases, for example, heritability is typically concentrated around immunological regions, like  the MHC/HLA region.

The issue of the implausibility of isotropic pleiotropy  is partially resolved by stratified linkage disequilibrium score regression[@finucane2015partitioning] (S-LDSC), an extension proposed by the same authors who devised LDSC.  S-LDSC allows the use of a pre-specified functional partitioning of the genome. While heritability is still assumed to be evenly distributed within a given partition, S-LDSC allows it to differ across partitions.


[//]: # (So it would be fair to say that while LDSC's key assumption is valid in a rough sense, it is not accurate in a granular sense.)

[//]: # (To discuss: Plietotopy an important fact from GWAS.  But not all variants equally likely.  e.g. HLA-> autoimmune disease)


## Usage
LDSC in its original form has two main usages

1. Estimation of the parameter $h^2$, the [heritability](Heritability.md) of phenotype of interest.  An estimate of this parameter can be directly read from the LDSC regression.
2. Detection of population stratification.



[//]: # ()
[//]: # ([//]: # Notes: Pleiotropy is key assumption.  Basically heritability is distributed more-or-less evening across the genome.  Contrast: MHC for inflammatory diseases.  Monogenic diseases.)

\bibliography

[//]: # (Main Reference:)

[//]: # (Bulik-Sullivan, Brendan K., et al. "LD Score regression distinguishes confounding from polygenicity in genome-wide association studies." Nature genetics 47.3 &#40;2015&#41;: 291-295.)