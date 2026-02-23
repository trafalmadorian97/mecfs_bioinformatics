---
hide:
- navigation
- toc
---
# CT-LDSC
Cross Trait Linkage Disequilibrium Score Regression (CT-LDSC)[@bulik2015atlas] is an extension of [Linkage Disequilibrium Score Regression](LDSC.md) (LDSC)[@bulik2015ld] that estimates [genetic correlation](Genetic_Correlation.md).

## Derivation

### Data Generating Model

Recall that LDSC assumes a linear data-generating model. CT-LDSC assumes two linear data-generating models, one for each of two traits:

$$
\begin{align}
y_1&= Y\beta + \delta \label{dg1}\\
y_2&= Z\gamma + \epsilon \label{dg2}
\end{align}
$$


where:

- There are $M  \gg 0$ genetic variants.  
- There are $N_1 \gg 0$ individuals in the GWAS of trait 1 and $N_2 \gg 0$ individuals in the GWAS of trait 2.
- There are $N_s$ individuals included in both GWAS.  Without loss of generality, assume these $N_s$ individuals are listed first in the lists of participants in both studies.
- $y_1\in \mathbb{R}^{N_1}$ and $y_2 \in \mathbb{R}^{N_2}$ are the vectors of phenotypes for the GWAS of trait 1 and the GWAS of trait 2, respectively.
- $Y\in\mathbb{R}^{N_1\times M}$ and $Z\in\mathbb{R}^{N_2\times M}$ are the genotype matrices from the two GWAS, normalized to have columns with sample mean 0 and variance 1.  $Y_j\in\mathbb{R}^{N_1}$ and $Z_j \in \mathbb{R}^{N_2}$ denote the $j$th columns of the two matrices.
- $\beta,\gamma\in\mathbb{R}^M$ are the vectors of true genetic effect sizes of each genetic variant on the two traits. 
- $Y\beta\in\mathbb{R}^{N_1}$ and $Z\gamma\in\mathbb{R}^{N_2}$ are thus the vectors of genetic effects in the two GWAS.
- $\epsilon\in\mathbb{R}^{N_1}$ and $\delta\in\mathbb{R}^{N_2}$ are the vectors of non-genetic effects in the two GWAS.

We model $Y,Z,\beta,\gamma, \delta,\epsilon$ as random variables with the following properties:


$$
\begin{align}
\mathbb{Var(\delta)} &= (1-h_1^2)I & \text{ for some }h_1>0\\
\mathbb{Var(\epsilon)} &= (1-h_2^2)I& \text{ for some }h_2>0\\
\mathbb{Cov}(\delta_j,\epsilon_k) &= \begin{cases}
\rho_e &\text{ if } j=k \le N_s \\
0 & \text{ else }
\end{cases} & \text{ for some }\rho_e>0 \label{cov_delta_epsilon}\\
\mathbb{Var}(\beta) &= \frac{1}{M} h_1^2 I \label{var_beta} \\
\mathbb{Var}(\gamma) &= \frac{1}{M} h_2^2 I  \label{var_gamma}  \\
\mathbb{Cov}(\beta, \gamma)&= \frac{1}{M}\rho_gI & \text{ for some }\rho_g>0 \label{cov_beta_gamma}\\
\mathbb{E}(\delta)&=0\label{e_delta}\\
\mathbb{E}(\beta)&=0  \label{e_beta} \\
\mathbb{E} Y &=0 \label{e_Y} \\ 
\mathbb{E}(\epsilon)&=0   \label{e_epsilon} \\
\mathbb{E}(\gamma)&=0  \label{e_gamma} \\
\mathbb{E} (Z) &=0  \label{e_Z} \\ 
\mathbb{E}(Y_i^2) &=\mathbb{E}(Z_i^2)=1 &\text{ for all }i \label{y_z_var}
\end{align}
$$

- Assumption ($\ref{cov_delta_epsilon}$) specifies that error is uncorrelated between GWAS, except for individuals who are participants both studies.
- Assumptions ($\ref{var_beta}$), ($\ref{var_gamma}$), and ($\ref{cov_beta_gamma}$) specify an isotropically polygenic architecture, analogous to the one used in [LDSC](LDSC.md).
- Assumptions ($\ref{e_delta}$), ($\ref{e_beta}$), ($\ref{e_Y}$), ($\ref{e_Z}$), and ($\ref{y_z_var}$) specify that the data has been pre-normalized.




Furthermore, we assume the following relationships between the random variables, which are faily standard in linear data-generating models:

- The rows of $Y$ and $Z$ are mutually independent, except in the case that a row of $Y$ and $Z$ refer to the same individual (one of the $N_s$ individuals present in both GWAS).
- The rows of $Y$ and $Z$ are identically distributed.
- $Y, \beta, \delta$ are all mutually independent, as are $Z,\gamma, \epsilon$.
- $\beta$ is independent of $\epsilon$ and $Z$, and $\gamma$ is independent of $\delta$ and $Y$.


Note that this implies that the true genetic effects have mean zero:

$$
\begin{align}
\mathbb{E}( Y\beta) &= \mathbb{E}(Y) \mathbb{E}(\beta) &\text{ by independence}\\
&=0
\end{align}
$$

and similarly, $\mathbb{E}(Z\gamma)=0$.

Define the following quantities related to Linkage Disequilibrium (LD):

- The LD between SNP $j$ and SNP $k$ is denoted by $r_{jk}:=\mathbb{E}(Y_{i,j}Y_{i,k})=\mathbb{E}(Z_{q,j}Z_{qk})$, (which does not depend on the individuals $i$ and $q$ by our assumption that the rows of $Y$ and $Z$ are iid).

[//]: # (- The empirical LD between $j$ and $k$ is defined to be $\tilde{r}_{jk}:=\frac{1}{N}X_{:,j}^T X_{:,k}$.)
- The LD score of a SNP $j$ is defined to be $l_j:= \sum_k r_{jk}^2$. It quantifies the magnitude of the dependence between $j$ and other SNPs.


### Genetic Covariance

First, let us compute the genetic covariances between the two phenotypes.  Let $X\in\mathbb{R}^M$ denote the genotype of an arbitrary individual.  By definition, genetic covariance is:

$$
\begin{align}
&\mathbb{Cov}\left( \sum_j X_j \beta_j, \sum_k X_k \gamma_k \right)\\
&= \mathbb{E}\left(  \left(\sum_j X_j \beta_j\right)\left(\sum_k X_k \gamma_k \right)  \right)& \text{mean 0}\\
&=\sum_{j,k} \mathbb{E} \left( X_j X_k \beta_j \gamma_k \right) & \text{linearity}\\
&=\sum_{j,k}\mathbb{E} (X_jX_k) \mathbb{E}(\beta_j \gamma_k) & \text{independence}\\
&=\sum_{j}\mathbb{E} (X_j^2) \mathbb{E}(\beta_j \gamma_j) & \text{by (\ref{cov_beta_gamma})}\\
&=\sum_j \mathbb{E} (\beta_j \gamma_j) & \text{ by (\ref{y_z_var})}\\
&=\rho_g. & \text{by (\ref{cov_beta_gamma})}
\end{align}
$$

### Genetic Correlation

Note that the model ($\ref{dg1},\ref{dg2}$) is an extension of the model used in [LDSC](LDSC.md).  By the derivation of LDSC, we have $\mathbb{Var}(\sum_j X_j \beta_j)=h_1^2$ and $\mathbb{Var}(\sum_j X_j \gamma_j)=h_2^2$.

The genetic correlation of the two traits can be computed as their genetic covariance divided by the square root of the product of their heritabilities:

$$
\frac{\rho_g}{h_1h_2}.
$$

### Regression equation

Using the same logic as in [derivation of LDSC](LDSC.md), we approximate the Wald $\chi^2$ statistics of a variant $j$ in the two GWAS: 

$$
\begin{align}
\chi_{j,1}^2 &\approx N \hat{\beta}_{j,1}^2=\frac{(y_1^T Y_{j})^2}{N_1}\\
\chi_{j,2}^2 &\approx N \hat{\beta}_{j21}^2=\frac{(y_1^T Z_{j})^2}{N_2}.
\end{align}
$$

It follows that the corresponding z-statistics are: 

$$
\begin{align}
z_{j,1} &\approx \frac{(y_1^T Y_{j})}{\sqrt{N_1}} \label{d1}\\
z_{j,2} &\approx \frac{y_1^T Z_{j}}{\sqrt{N_2}} \label{d2}.
\end{align}
$$


Recall that in [LDSC](LDSC.md), the regression dependent variable was the Wald $\chi^2$ statistic.  In CT-LDSC, the regression dependent variable is the product of the $z$-statistics.  To derive the regression equation, we must estimate the expectation of this product: $\mathbb{E}(z_{j,1}z_{j,2})$.  By the Tower Law of Expectation (see Grimmett and Stirzaker[@grimmett2020probability]), 

$$
\begin{align}
\mathbb{E}(z_{j,1}z_{j,2})= \mathbb{E}\left( \mathbb{E}(z_{j,1}z_{j,2}|Z,Y) \right). \\
\end{align}
$$  


Evaluating the inner expectation,

$$
\begin{align}
&\mathbb{E}(z_{j,1}z_{j,2}|Z,Y)\\
&\approx\frac{1}{\sqrt{N_1N_2}} \mathbb{E}\left(Y_j^T y_1 y_1^T Z_j|Z,Y \right) & \text{ by (\ref{d1}, \ref{d2})}\\
&=\frac{1}{\sqrt{N_1N_2}}Y^T_j \mathbb{E}\left( (Y\beta+\delta)(Z\gamma+\epsilon)^T  |Z,Y\right)Z_j & \text{ by (\ref{dg1},\ref{dg2})}\\
&=\frac{1}{\sqrt{N_1N_2}}Y^T_j \left(Y \mathbb{E}(\beta \gamma^T)Z^T +  \mathbb{E}(\delta \gamma^T) Z^T + Y \mathbb{E}(\beta \epsilon^T) + \mathbb{E}(\epsilon \delta^T) \right)Z_j & \text{linearity, independence}\\
&=\frac{1}{\sqrt{N_1N_2}}Y^T_j\left(\frac{\rho_g}{M}YZ^T + \rho_e Q\right)Z_j & \text{ by (\ref{cov_beta_gamma},\ref{cov_delta_epsilon})}\\
&=\frac{1}{\sqrt{N_1N_2}}\left(\frac{\rho_g}{M}Y_j^TYZ^TZ_j + \rho_e Y_j^T QZ_j\right) \label{cond_exp}
\end{align}
$$

where $Q\in\mathbb{R}^{N_1\times N_2}$ and 

$$
\begin{align}
Q_{i,j}&= \begin{cases}
1 & \text{ if }i=j \le N_s\\
0 & \text{else}
\end{cases}
\end{align}
$$

Our goal is to derive an expression for the unconditional expectation of $z_{j,1}z_{j,2}$.  To achieve this, we take the expectation of both terms in $(\ref{cond_exp})$.  First, we have

$$
\begin{align}
&\mathbb{E}(Y_j^TYZ^TZ_j)\\
&= Y_j^T \left( \sum_{k=1}^M Y_{i,k} Z_{q,k} \right)_{(i,q)}  Z_j \\
&=\mathbb{E} \sum_{i=1}^{N_1}\sum_{q=1}^{N_2}\sum_{k=1}^M\left( Y_{i,j}Y_{i,k} Z_{q,j} Z_{q,k}  \right) & \text{ def of matrix product}
\end{align}
$$

There are $N_1N_2$ possible pairs of the indexes $i$ and $q$.  For each such pair, there are two cases:

- Case 1: $i$ and $q$ refer to the same individual.  There are $N_s$ such index pairs. In this case we have
 
$$
\begin{align}
&\sum_{k=1}^M \mathbb{E}(Y_{i,j} Y_{i,k} Z_{q,j} Z_{q,k} )\\
&=\sum_{k=1}^M \mathbb{E}(Y_{i,j}^2 Y_{i,k}^2)\\
&\approx  \sum_{k=1}^M (1 + 2r_{j,k}^2)  \label{isserlis_approx_line}\\
&=M + 2\sum_{k=1}^M r_{j,k}^2\\
&= M + 2l_j 
\end{align}
$$

Where we have approximated the random variables as having a normal distribution, and used 
[ Isserlis's Theorem](https://en.wikipedia.org/wiki/Isserlis%27s_theorem).  Recall that we also used Isserlis's theorem at a similar point in the original derivation of [LDSC](LDSC.md#expectation-of-empirical-ld-scores).


- Case 2: $i$ and $q$ refer to different individuals.  There are $N_1N_2-N_s$ such pairs.

$$
\begin{align}
&\sum_{k=1}^M \mathbb{E}(Y_{i,j} Y_{i,k} Z_{q,j} Z_{q,k} )\\
&=\sum_{k=1}^M \mathbb{E}(Y_{i,j} Y_{i,k}) \mathbb{E}(Z_{q,j} Z_{q,k} ) & &\text{ independence}  \\
&=\sum_{k=1}^M r_{jk}^2\\
&= l_{j}
\end{align}
$$

Combining the two cases yields

$$
\begin{align}
&\mathbb{E}(Y_j^TYZ^TZ_j)\\
&=N_s(M+2l_j) + (N_1N_2-N_s)l_j\\
&=MN_s + (N_1N_2+N_s)l_j\\
&\approx MN_s + N_1N_2 l_j. \label{exp_prod_1}
\end{align}
$$


Returning to the second term in $(\ref{cond_exp})$, we have that 

$$
\begin{align}
\mathbb{E}(Y_j^TQZ_j)=N_s \label{exp_prod_2}
\end{align}
$$ 

by the assumption of independence of the genotypes of distinct individuals.

Substituting $(\ref{exp_prod_1})$ and $(\ref{exp_prod_2})$ into $(\ref{cond_exp})$ yields:

$$
\begin{align}
&\mathbb{E}(z_{j,1}z_{j,2})\\
&\approx \frac{1}{\sqrt{N_1N_2}}\left( \frac{\rho_g}{M}  (MN_s + N_1N_2 l_j) +\rho_e N_s \right)\\
&= \frac{\rho_g\sqrt{N_1N_2}}{M}l_j + \frac{N_s}{\sqrt{N_1N_2}    } (\rho_g + \rho_e   ).\label{ct_ldsc_eqn}
\end{align}
$$

$(\ref{ct_ldsc_eqn})$ is the key regression equation in CT-LDSC.


### Note on derivation

The derivation in **the supplementary material to the original paper[@bulik2015atlas]** contains several typos and implicit approximations.  In the version above, I have corrected typos and make approximations explicit.