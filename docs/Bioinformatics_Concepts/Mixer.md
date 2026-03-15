# MiXeR

MiXeR[@holland2020beyond;@frei2019bivariateMkd] is a family of parametric probabilistic models for the distribution of variant effect sizes in one or more GWAS.  Authors will often fit MiXeR models to summary statistics from traits of interest, then inspect the resulting model parameters to understand the genetic architecture of these traits.

The authors of MiXeR have published numerous papers proposing variants of the core MiXeR model, each of which is suitable for slightly different circumstances.

##  Univariate MiXeR

Like [LDSC](LDSC.md), MiXeR assumes a linear data generating model:


$$
y = \sum_{i}\beta_i g_i +\epsilon
$$

where 


- $y\in\mathbb{R}^N$ is the vector of phenotypes of study participants.
- $\beta_i\in\mathbb{R}$ is the true causal regression coefficient of the $i$th genetic variant.
- $g_i \in\mathbb{R}^N$ is the vector of genotypes of study participants at the $i$th genetic variant.  
- $\epsilon\in\mathbb{R}^N$ is the vector of environmental and other effects.

MiXeR models the true causal effect sizes via a mixture distribution.  For any variant $i$, we have

$$
\beta_i \sim \pi_1 \mathcal{N}(0, \sigma^2_\beta) + (1-\pi_1) \mathcal{N}(0,0)
$$

Where:

- $\pi_1\in (0,1)$ is the proportion of genetic variants that are causal for the trait of interest.  
- $\mathcal{N}(0,\sigma_\beta^2)$ is a Gaussian distribution with variance $\sigma_\beta^2>0$, and $\mathcal{N}(0,0)$ is a Dirac delta distribution[^dirac_note].

Moreover, we assume that for genetic variants $i\ne j$, $\beta_i$ and $\beta_j$ are independent.


[^dirac_note]: $\beta_i\sim \mathcal{N}(0,0)$ just means that $P(\beta=0)=1$.