# Colocalization

Colocalization is a biostatistical technique to evaluate the evidence that the GWAS signal at a particular locus for two different traits is driven by the same underlying causal variant.  

Because it evaluates commonality between GWAS signals, colocalization is somewhat similar in spirit to [genetic correlation](Genetic_Correlation.md) techniques like [LDSC](LDSC.md).  However, whereas genetic correlation evaluates global polygenic overlap, colocalization is narrowly focused on a single locus.  It is thus possible for two GWAS to strongly colocalize at one locus, but not at another. Colocalization techniques enable a nuanced understanding of where two traits share biological pathways and where they do not.


Colocalization is often paired with [Mendelian Randomization](Mendelian_Randomization.md), another local, variant-focused method.


## Bayes factor colocalization

Bayes factor colocalization[@giambartolomei2014bayesian] is elegant and simple, but requires certain restrictive assumptions.  It is used as a basis for more advanced colocalization methods.



### Data generating model

We start by assuming we are working with quantitative traits.  Bayes factor colocalization assumes the following generative model.


- There are $Q$ genetic variants at the locus of interest.
- We focus on two GWAS.  There are $N_1$ participants in the first GWAS and $N_2$ participants in the second GWAS.
- $y_k\in\mathbb{R}^{N_k}$ is the vector of phenotypes for participants in GWAS $k$.
- $G_k \in \mathbb{R}^{N_k \times Q}$ is the genotype matrix for GWAS $k$ at the locus of interest.
- $g_{k,i}$ is the $i$th column of $G_k$.
- $\beta_{k,i}$: The causal, joint underlying regression coefficient of variant $i$ for trait $k$.
- $\hat\beta_{k,i}$: The GWAS marginal regression coefficient of variant $i$ in GWAS $k$.
- $V_{k,i}$: The squared standard error of $\hat\beta_{k,i}$.
- $\epsilon_k\in\mathbb{R}^{N_k}$:  The part of trait $k$ that is not determined by variants at the locus of interest.  Thus $\epsilon_k$ includes both environmental effects and the effects of variants at other loci.
- $c_k$ is the number of covariates used in the trait $k$ model.
- $X_k\in\mathbb{R}^{N_k\times c_k}$ is a matrix of trait $k$ covariates.
- $\alpha_k\in\mathbb{R}^{c_k}$ is the vector of true trait $k$ covariate coefficients.


Our data generating model is:

$$
\begin{align}
y_k &= X_k\alpha_k +G_k \beta_k + \epsilon_k\\
\epsilon_k &\sim N(0,\sigma^2_k I).
\end{align}
$$



[//]: # (Furthermore, we assume that there is at most one causal variant at the locus for each of the two traits.  Let $S_1,S_2 \in \{ 0,\ldots, Q\}$ be the integer index of the causal variant for each trait if there is one, and 0 if there is no causal variant.  Thus)

[//]: # ()
[//]: # ($$)

[//]: # (\begin{align})

[//]: # (\beta_{k,i} &=0     &\text{ if }i\ne S_k,\\)

[//]: # (\beta_{k,i}&\ne 0  &\text{ if }i=S_k.)

[//]: # (\end{align})

[//]: # ($$)


We assume:

- $\epsilon_1$ and $\epsilon_2$ are independent.  This assumption is usually satisfied by assuming that there are no overlapping participants between the two studies.

[//]: # (- We assume that $\beta_{k,i}| | \{ S_k =i \} \sim \mathcal{N}&#40;0,\omega_k^2&#41;$ .)


### Category vectors

A variant $i$ is said to be causal for trait $k$ if $\beta_{k,i}\ne 0$.  Use $\mathcal{C}=\{0,1,2,12\}$  to express the 4 disjoint possibilities for any given variant:

- $0$: The variant is not causal for either trait.
- $1$: The variant is causal for trait 1 only.
- $1$: The variant is causal for trait 2 only.
- $12$: The variant is causal for both trait 1 and trait 2.

Call $\mathcal{C}^Q$ the set of category vectors.  Each category vector corresponds to a decision about the causal nature of each variant at the locus.  Assume a prior under which the causality of the $Q$ variants are independent, and the probabilities are of the four possibilities for any given variant are $p_0,p_1,p_2,p_{1,2}$.

Thus the prior probability of any category vector $c$ is 

$$
\begin{align}
P(c) = \prod_{i=1}^Q p_{c_i}.
\end{align}
$$


### The single causal variant assumption


Let $\mathcal{A}$ denote the event that at most one variant is causal for each of the two traits.  Define the function $\overline{S}$, whose domain is the set of category vectors corresponding to $\mathcal{A}$.  Define $\overline{S}(c):=(S_1,S_2)$, where $S_1$ is the index of the causal variant for trait 1, or 0 if there is no causal variant, and $S_2$ is the index of the causal variant for trait $2$, or 0 if there is no causal variant. A pair $(S_1,S_2)$ is called a configuration. There are $(Q+1)^2$ configurations.
 
Bayes factor colocalization operates conditional on  $\mathcal{A}$, the assumption that at most one variant is causal.  This assumption is unrealistic for some loci, but permits significant simplification.

[//]: # (For an arbitrary variant $i$, define)

[//]: # ()
[//]: # (- $p_0$: the probability that variant $i$ is not causal for either of the two traits.)

[//]: # (- $p_1$: the probability that variant $i$ is causal for trait 1 only.)

[//]: # (- $p_2$: the probability that variant $i$ is causal for trait 2 only.)

[//]: # (- $p_{1,2}$: the probability that variant $i$ is causal for both trait 1 and trait 2.)

[//]: # ()
[//]: # ()
[//]: # (Define  $\mathcal{C}$ denote the set )

[//]: # ()

## Hypotheses

[//]: # (HypothesesA realization of $&#40;S_1,S_2&#41;$ defining the causal variant for each of the two traits is called a configuration. Since there are $Q+1$ possible values of each of the $S_k$, there are $&#40;Q+1&#41;^2$ configurations.  )

We define 5 hypothesis, which partition the configurations, as follows:

- $H_0$: Corresponds to the configurations $\mathcal{S}_0=\{ (0,0) \}$.  No causal variants at the locus.
- $H_1$: Corresponds to the configurations $\mathcal{S}_1=\{ (i,0) : i \ge  1 \}$.  A causal variant for trait 1 but not trait 2 at the locus.
- $H_2$: Corresponds to the configurations $\mathcal{S}_2=\{ (0,j) : j \ge  1 \}$.  A causal variant for trait 2 but not trait 1 at the locus.
- $H_3$: Corresponds to the configurations $\mathcal{S}_3=\{ (i,j) : i,j \ge  1 ,i\ne j\}$. Both traits have causal variants at the locus, but they are different .
- $H_4$: Corresponds to the configurations $\mathcal{S}_4=\{ (i,i) : i \ge  1\}$. Both traits have causal variants at the locus, but they are different .

See the illustration below,

### Configuration odds


$$
\begin{align}
\rho_1&:=\frac{p_1}{p_0}\\
\rho_2&:= \frac{p_2}{p_0}\\
\rho_{12}&:=\frac{p_{12}}{p_0}.
\end{align}
$$

For some $h\in \{0,1,2,3,4\}$, let $S\in\mathcal{S}_h$. By the prior, we have


$$
\begin{align}
P(S)&=
\begin{cases}
p_0^Q & h=0\\
p_1 p_{0}^{Q-1} & h=1\\
p_2 p_{0}^{Q-1} & h=2\\
p_1 p_2 p_{0}^{Q-2} & h=3\\
p_{12} p_{0}^{Q-1} & h=2\\
\end{cases} \label{unconditional_prob}
\end{align}
$$



Let $S^0:=(0,0)$. We have


$$
\begin{align}
&\frac{P(S|\mathcal{A})}{P(S^0|\mathcal{A})}\\
&=\frac{P(S \cap \mathcal{A}) / P(\mathcal{A})}{P(S^0 \cap \mathcal{A})/ P(\mathcal{A})}\\
&=\frac{P(S)}{P(S^0)}\\
&=\frac{P(S)}{p_0^Q}\\
&=\begin{cases}
1 & h=0\\
\rho_1 & h=1\\
\rho_2 & h=2\\
\rho_1 \rho_2 & h=3\\
\rho_{12} & h=4
\end{cases} & \text{by (\ref{unconditional_prob})}
\end{align}
$$



## SUSIE colocalization

todo


## Links


[Talk on colocalization at International Statistical Genetics Workshop](https://youtu.be/QNUdAFn1Vas?si=dW0v2QJF0kzZIIea)