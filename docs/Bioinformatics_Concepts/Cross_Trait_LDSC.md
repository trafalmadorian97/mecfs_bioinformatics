# Cross Trait Linkage Disequilibrium Score Regression
Cross Trait Linkage Disequilibrium Score Regression (CT-LDSC)[@bulik2015atlas] is an [extension of Linkage Disequilibrium Score Regression](LDSC.md) (LDSC) that allows the estimation of [genetic correlation](Genetic_Correlation.md).

## Derivation

### Data Generating Model

Recall that LDSC assumes a linear data-generating model for the trait of interest.  Similarly, CT-LDSC assumes two linear data-generating models.

$$
\begin{align}
y_1&= Y\beta + \delta\\
y_2&= Z\gamma + \epsilon
\end{align}
$$


where:

- There are $M  \gg 0$ genetic variants.  
- There are $N_1 \gg 0$ individuals in the GWAS of trait 1 and $N_2 \gg 0$ individuals in the GWAS of trait 2.
- There are $N_s$ individuals included in both GWAS.  These $N_s$ individuals are listed first in the list of participants in both studies.
- $y_1\in \mathbb{R}^{N_1}$ and $y_2 \in \mathbb{R}^{N_2}$ are the vectors of phenotypes for trait 1 and trait 2 respectively.
- $Y\in\mathbb{R}^{N_1\times M}$ and $Z\in\mathbb{R}^{N_2\times M}$ are the genotype matrices from the two GWAS, normalizd to have columns with sample mean 0 and variance 1.
- $\beta,\gamma\in\mathbb{R}^M$ are the vectors of true genetic effect sizes for each genetic variant for the two traits. 
- $Y\beta\in\mathbb{R}^{N_1}$ and $Z\gamma\in\mathbb{R}^{N_2}$ are thus the vectors of true genetic effects for all individual in the two GWAS.
- $\epsilon\in\mathbb{R}^{N_1} and $\delta\in\mathbb{R}^{N_2}$ are the vectors of non-genetic effects for all individuals in the two GWAS.

Furthermore, we model $Y,Z,\beta,\gamma, \delta,\epsilon$ as random variables with the following properties


$$
\begin{align}
\mathbb{Var(\delta)} &= (1-h_1^2)I\\
\mathbb{Var(\epsilon)} &= (1-h_2^2)I\\
\mathbb{Cov}(\delta_j,\epsilon_k) &= \begin{cases}
\rho_e &\text{ if } j=k \le N_s \\
0 & \text{ else }
\end{cases}\\
\mathbb{Var}(\beta) &= \frac{1}{M} h_1^2 I\\
\mathbb{Var}(\gamma) &= \frac{1}{M} h_2^2 I\\
\mathbb{Cov}(\beta, \gamma)&= \frac{1}{M}\rho_gI \label{cov_beta_gamma}\\
\mathbb{E}(\delta)&=0\\
\mathbb{E}(\beta)&=0 \\
\mathbb{E} Y &=0 \\ 
\mathbb{E}(\epsilon)&=0\\
\mathbb{E}(\gamma)&=0 \\
\mathbb{E} (Z) &=0 \\ 
\mathbb{E}(Y_i^2) &=\mathbb{E}(Z_i^2)=1 &\text{ for all }i \label{y_z_var}
\end{align}
$$

- The rows of $Y$ and $Z$ are mutually independent, except in the case that a row of $Y$ and $Z$ refer to the same individual (one of the $N_s$ individuals present in both GWAS).
- Furthermore, $Y, \beta, \delta$ are all mutually independent, as are $Z,\gamma, \epsilon$
- $\beta$ is independent of $\epsilon$ and $Z$, and $\gamma$ is independent of $\delta$ and $Y$


Note that this implies that the true genetic effects have mean zero:

$$
\begin{align}
\mathbb{E}( Y\beta) &= \mathbb{E}(Y) \mathbb{E}(\beta) &\text{ by independence}\\
&=0
)
\end{align}
$$

and similarly, $\mathbb{E}(Z\gamma)=0$

Furthermore, define the following quantities related to Linkage Disequilibrium (LD).

- The LD between SNP $j$ and SNP $k$ is denoted by $r_{jk}:=\mathbb{E}X_{i,j}X_{i,k}$, (which does not depend on the individual $i$ by our assumption that the rows of $X$ are iid).
- The empirical LD between $j$ and $k$ is defined to be $\tilde{r}_{jk}:=\frac{1}{N}X_{:,j}^T X_{:,k}$.
- The LD score for a SNP $j$ is defined to be $l_j:= \sum_k r_{jk}^2$.


### Genetic Covariance

First, let us compute the genetic covariances between the two phenotypes.  Let $X\in\mathbb{R}^M$ denote the genotype of an arbitrary individual.  

$$
\begin{align}
&\mathbb{Cov}\left( \sum_j X_j \beta_j, \sum_k X_k \gamma_k \right)\\
&=\sum_{j,k} \mathbb{E} \left( X_j X_k \beta_j \gamma_k \right) & \text{mean 0}\\
&=\sum_{j,k}\mathbb{E} (X_jX_k) \mathbb{E}(\beta_j \gamma_k) & \text{independence}\\
&=\sum_{j}\mathbb{E} (X_j^2) \mathbb{E}(\beta_j \gamma_j) & \text{by (\ref{cov_beta_gamma})}\\
&=\sum_j \mathbb{E} (\beta_j \gamma_j) & \text{ by (\ref{y_z_var})}\\
&=\rho_g & \text{by (\ref{cov_beta_gamma})}
\end{align}
$$

### Genetic Correlation

Note that this model includes the model used in the derivation of [LDSC](LDSC.md).  Thus by the logic of that derivation, we have $\mathbb{Var}(\sum_j X_j \beta_j)=h_1^2$ and $\mathbb{Var}(\sum_j X_j \gamma_j)=h_2^2$.

The genetic correlation of the two traits is

$$
\frac{\rho_g}{h_1h_2}
$$

### Regression equation

Using the same logic as in [derivation of LDSC](LDSC.md), we can derive approximations for the Wald $\chi^2$ statistics for the GWAS regressions of traits 1 and 2 on variant $j$: 

$$
\begin{align}
\chi_{j,1}^2 &\approx N \hat{\beta_{j,1}^2}=\frac{(y_1^T Y_{j})^2}{N}\\
\chi_{j,2}^2 &\approx N \hat{\beta_{j21}^2}=\frac{(y_1^T Z_{j})^2}{N}.
\end{align}
$$

From this, it follows that the z statistics for the GWAS regressions on variant $j$ are approximately

$$
\begin{align}
\z_{j,1} &\approx \frac{(y_1^T Y_{j})}{\sqrt{N}}\\
\z_{j,2}^2 &\approx \frac{y_1^T Z_{j}}{\sqrt{N}}.
\end{align}
$$


