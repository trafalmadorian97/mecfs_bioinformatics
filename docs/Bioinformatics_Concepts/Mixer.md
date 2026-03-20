---
hide:
- toc
---
# MiXeR

MiXeR[@holland2020beyond;@frei2019bivariate] is a parametric model for the distribution of GWAS effect sizes.  The parameters of fitted MiXeR models provide insights into genetic architecture.  The three main versions of MiXeR are univariate MiXeR[@holland2020beyond], bivariate MiXeR[@frei2019bivariate], and GSEA MiXeR[@frei2022improved].

[//]: # (The authors of MiXeR have published numerous papers proposing variants of the core MiXeR model, each of which is suitable for slightly different circumstances.)

##  Univariate MiXeR

Univariate MiXeR[@holland2020beyond] models effect sizes in a GWAS of a single trait.

Like [LDSC](LDSC.md), MiXeR assumes a linear data generating model:


$$
\begin{align}
y = \sum_{i}\beta_i g_i +\epsilon \label{mixer_dgm}
\end{align}
$$

where 


- $y\in\mathbb{R}^N$ is the vector of phenotypes of study participants.  $y$ is normalized to have sample mean 0 and sample variance 1.
- $\beta_i\in\mathbb{R}$ is the true causal regression coefficient of the $i$th genetic variant.
- $g_i \in\mathbb{R}^N$ is the vector of genotypes of study participants at the $i$th genetic variant.  $g_i$ is normalized to have sample mean 0.
- $\epsilon\in\mathbb{R}^N$ is the vector of environmental effects.


Similar to LDSC, we view all of these variables as random, and consider the $N$ study participants to be independent and identically distributed. We have:

$$
\begin{align}
\mathrm{E}(y_j)&=0 & \text{ for all $j$}\\
\mathrm{Var}(y_j)&=1& \text{ for all $j$}\\
\sum_i \mathrm{Var}(\beta_i) &= h^2 &\text{ $h^2$ is heritability}\\
\mathrm{E} \epsilon &= 0\\
\mathrm{Var} (\epsilon) &= 1 -h^2\\
\end{align}
$$


Furthermore,
 
- Define the population variance of genetic variant $i$:

$$
\begin{align}
H_i&:= \mathrm{Var}(g_{j,i})
\end{align}
$$



for any individual $j$, 

- Define the sample variance of genetic variant $i$,

$$
\begin{align}
\hat H_i := \frac{\lVert g_i \rVert^2}{N}. \label{mixer_h_def}
\end{align}
$$

.

- Define the population correlation between variants $i$ and $j$

$$
\begin{align}
r_{i,j}:=\mathrm{Corr}(g_{k,i},g_{k,j})
\end{align}
$$

for any individual $k$.




- Define 

$$
\begin{align}
\hat r_{i,j}&:= \frac{g_i^T g_j}{\lVert g_i \rVert \lVert g_j \rVert  }\\
&=\frac{g_i^Tg_j}{N \sqrt{\hat H_i \hat H_j} } \label{corr_def},
\end{align}
$$

the sample correlation of variants $i$ and $j$.


MiXeR assumes true causal effect sizes follow a mixture distribution.  For any variant $i$, we have

$$
\begin{align}
\beta_i \sim C \mathcal{N}(0, \sigma^2_\beta) + (1-C) \mathcal{N}(0,0) \label{uni_mixer_core}
\end{align}
$$

Where:


- $C$ is a Bernoulli random variable with parameter $\pi_1$, where $\pi_1\in (0,1)$ is the proportion of genetic variants that are causal for the trait of interest.  
- $\mathcal{N}(0,\sigma_\beta^2)$ is a Gaussian distribution with variance $\sigma_\beta^2>0$, and $\mathcal{N}(0,0)$ is a Dirac delta distribution[^dirac_note].
- The Bernoulli and Gaussian distributions are independent of one another.

Moreover, for genetic variants $i\ne j$, $\beta_i$ and $\beta_j$ are independent.

By allowing some variants to have no effect, $(\ref{uni_mixer_core})$ is more general than the commonly used infinitesimal genomic model, which assumes that all variants affect the phenotype.


###  Distribution of $z$-scores

$(\ref{mixer_dgm})$ and $(\ref{uni_mixer_core})$ specify the core univariate MiXeR model. The next steps are:

- Write this model as a probability distribution over $z$-scores of genetic variants.  
- Derive an efficient algorithm to fit the MiXeR model by maximum likelihood, given GWAS $z$-scores.

We begin by deriving an expression for $\hat{\beta}_i$, the univariate regression coefficient of variant $i$.


$$
\begin{align}
\hat{\beta}_i&=\frac{y^Tg_i}{g_i^T g_i} \text{ (OLS coefficient formula)}\\
&= \frac{y^Tg_i}{N \hat H_i} \text{ (By (\ref{mixer_h_def}))}\\
&=\frac{(\sum_j \beta_j g_j + \epsilon)^Tg_i}{N \hat H_i} \text{ (By (\ref{mixer_dgm}))}\\
&=\sum_j \beta_j\frac{g_j^T g_i}{N \hat H_i }  + \frac{\epsilon^T g_i}{N \hat H_i}\\
&=\sum_j \beta_j\frac{\hat r_{i,j}\sqrt{\hat H_i \hat H_j} N  }{N\hat H_i }  + \frac{\epsilon^T g_i}{N \hat H_i} \text{ (by (\ref{corr_def}))}\\
&=\sum_j \beta_j\hat r_{i,j}\sqrt{\frac{\hat H_j}{\hat H_i}} + \frac{\epsilon^T g_i}{N \hat H_i}
\end{align}
$$

Next, we compute the standard error of $\hat{\beta_i}$.

$$
\begin{align}
\mathrm{SE}(\hat\beta_i)&\approx \sqrt{\frac{\frac{1}{N}\lVert y- \hat{\beta}_i g_i \rVert^2}{g_i^T g_i}} & \text{ OLS SE formula}\\
&\approx \sqrt{\frac{\lVert y \rVert^2/N}{g_i^T g}} \label{se_step}\\
&= \sqrt{\frac{1}{g_i^T g}}& \text{ Normalization}\\
&=\frac{1}{\sqrt{N \hat H_i}} & \text{ By (\ref{mixer_h_def})}
\end{align}
$$

where $(\ref{se_step})$ follows since we assume we are studying a polygenic trait, for which no single variant can explain a significant proportion of the variance of the phenotype.


It follows that the $z$-score of variant $i$ is:

$$
\begin{align}
z_i&:=\frac{\hat{\beta}_i}{\mathrm{SE}(\hat\beta_i) }\\
&=\sum_j \beta_j \hat r_{i,j}\sqrt{\hat H_jN} + \frac{\epsilon^T g_i}{\sqrt{N \hat H_i}} \label{mixer_zscore_eq}
\end{align}
$$


$(\ref{mixer_zscore_eq})$ can be used to relate the distribution of $z$-scores (observed) to the distribution of the $\beta_i$ (implied by model parameters).

It can be useful to re-group the terms in $(\ref{mixer_zscore_eq})$. Define:

$$
\mathrm{LD}(i):=\{ j: r_{i,j=0} \}
$$

the set of variants not in linkage disequilibrium with $i$.


We have

$$
\begin{align}
z_i&=\underbrace{\sum_{j\in \mathrm{LD}(i)} \beta_j \hat r_{i,j}\sqrt{\hat H_jN}}_{=:\delta_i} +  \underbrace{\sum_{j\notin \mathrm{LD}(i)} \beta_j \hat r_{i,j}\sqrt{\hat H_jN} + \frac{\epsilon^T g_i}{\sqrt{N \hat H_i}}}_{=:e_i} 
\end{align}
$$

Note that even though the population correlation $r_{i,j}=0$ for all variants $j$ not  in $\mathrm{LD}(i)$, the sample correlation $\hat r_{i,j}$ will in general be nonzero.


### Variance decomposition





Next, let us compute the variance of $e_i$.


$$
\begin{align}
&\mathrm{Var}(e_i)\\
&=\mathrm{Var}\left( \sum_{j\notin \mathrm{LD}(i)}  \beta_j \hat r_{i,j}\sqrt{\hat H_jN} + \frac{\epsilon^T g_i}{\sqrt{N \hat H_i}}  \right) \\ 
&=
\end{align}
$$



to be continued

[^dirac_note]: $\beta_i\sim \mathcal{N}(0,0)$ means $P(\beta=0)=1$.