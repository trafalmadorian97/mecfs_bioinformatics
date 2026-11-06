# Colocalization

Colocalization is a biostatistical technique to evaluate the evidence that the GWAS signal at a particular locus for two different traits is driven by the same underlying causal variant.  

Because it evaluates commonality between GWAS signals, colocalization is somewhat similar in spirit to [genetic correlation](Genetic_Correlation.md) techniques like [LDSC](LDSC.md).  However, whereas genetic correlation evaluates global polygenic overlap, colocalization is narrowly focused on a single locus.  It is thus possible for two GWAS to strongly colocalize at one locus, but not at another. Colocalization techniques enable a nuanced understanding of where two traits share biological pathways and where they do not.


Colocalization is often paired with [Mendelian Randomization](Mendelian_Randomization.md), another local, variant-focused method.


## Bayes factor colocalization

Bayes factor colocalization is elegant and simple, but requires certain restrictive assumptions.  It is used as a basis for more advanced colocalization methods.



## Data generating model

We start by assuming we are working with quantitative traits.  Bayes factor colocalization assumes the following generative model.


- There are $Q$ genetic variants at the locus of interest.
- We focus on two GWAS.  There are $N_1$ participants in the first GWAS and $N_2$ participants in the second GWAS.
- $y_k\in\mathbb{R}^{N_k}$ is the vector of phenotypes for participants in GWAS $k$.
- $G_k \in \mathbb{R}^{N_k \times Q}$ is the genotype matrix for GWAS $k$ at the locus of interest.
- $g_{k,i}$ is the $i$th column of $G_k$.
- $\beta_{k,i}$: The true (joint) underlying coefficient of variant $i$ for the trait studied by GWAS $k$.
- $\hat\beta_{k,i}$: The GWAS marginal regression coefficient of variant $i$ in GWAS $k$.
- $V_{k,i}$: The squared standard error of variant $i$ in GWAS $k$.
- $\epsilon_k\in\mathbb{R}^{N_k}$:  The part of the phenotype for GWAS $k$ that is not determined by variants at the locus of interest.  Thus $\epsilon_k$ includes both environmental effects and the effects of variants at other loci.
- $c_k$ is the number of covariates used in GWAS $k$.
- $X_k\in\mathbb{R}^{N_k\times c_k}$ is a matrix of covariate values for GWAS $k$.
- $\alpha_k\in\mathbb{R}^{c_k}$ is the vector of covariate coefficients for GWAS $k$.


Our data generating model is:

$$
\begin{align}
y_k &= X_k\alpha_k +G_k \beta_k + \epsilon_k\\
\epsilon_k &\sim N(0,\sigma^2_k I).
\end{align}
$$



Furthermore, we assume that there is at most one causal variant at the locus for each of the two traits.  Let $S_1,S_2 \in \{ 0,\ldots, Q\}$ be the integer index of the causal variant for each trait if there is one, and 0 if there is no causal variant.  Thus

$$
\begin{align}
\beta_{k,i} &=0     &\text{ if }i\ne S_k,\\
\beta_{k,i}&\ne 0  &\text{ if }i=S_k.
\end{align}
$$


Additionally, we assume:

- $\epsilon_1$ and $\epsilon_2$ are independent conditional on $S_1,S_2$.  This assumption is usually satisfied by assuming that there are no overlapping participants between the two studies.
- We assume that $\beta_{k,i}| | \{ S_k =i \} \sim \mathcal{N}(0,\omega_k^2)$ .

Let us also define the per variant prior probabilities.  Fix an arbitrary variant $i$.  Define

- $p_0$: the probability that variant $i$ is not causal for either of the two traits.
- $p_1$: the probability that variant $i$ is causal for trait 1 only.
- $p_2$: the probability that variant $i$ is causal for trait 2 only.
- $p_{1,2}$: the probability that variant $i$ is causal for both trait 1 and trait 2.


## Configurations

A realization of $(S_1,S_2)$ defining the causal variant for each of the two traits is called a configuration.  Since there are $Q+1$ possible values of each of the $S_k$, there are $(Q+1)^2$ configurations.  We define 5 hypothesis, which partition the configurations, as follows:

- $H_0$: Corresponds to the configurations $\mathcal{S}_0=\{ (0,0) \}$.  No causal variants at the locus.
- $H_1$: Corresponds to the configurations $\mathcal{S}_1=\{ (i,0) : i \ge  1 \}$.  A causal variant for trait 1 but not trait 2 at the locus.
- $H_2$: Corresponds to the configurations $\mathcal{S}_2=\{ (0,j) : j \ge  1 \}$.  A causal variant for trait 2 but not trait 1 at the locus.
- $H_3$: Corresponds to the configurations $\mathcal{S}_3=\{ (i,j) : i,j \ge  1 ,i\ne j\}$. Both traits have causal variants at the locus, but they are different .
- $H_4$: Corresponds to the configurations $\mathcal{S}_4=\{ (i,i) : i \ge  1\}$. Both traits have causal variants at the locus, but they are different .

See the illustration below,

## Role of the prior



## SUSIE colocalization

todo


## Links


[Talk on colocalization at International Statistical Genetics Workshop](https://youtu.be/QNUdAFn1Vas?si=dW0v2QJF0kzZIIea)