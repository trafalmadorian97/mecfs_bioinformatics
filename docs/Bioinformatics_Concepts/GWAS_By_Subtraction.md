# GWAS by Subtraction

GWAS by subtraction[@demange2021investigating; @huang2024gwas] is a GenomicSEM[@grotzinger2019genomic] technique that orthogonally decomposes GWAS traits.

Here, we explain GWAS by subtraction twice: once at a high level via linear algebra, and more granularly via statistical modeling.

## Linear algebra

### Euclidian space

It is useful to understand GWAS-by-subtraction via linear algebra.


Consider a [Euclidian space](https://en.wikipedia.org/wiki/Euclidean_space) in which:

- GWAS traits are vectors.
- The [inner product](https://en.wikipedia.org/wiki/Inner_product_space) of two traits is their [genetic covariance](Genetic_Correlation.md).  Denote the inner product of $u$ and $v$ as $\langle u,v \rangle$.
- We assume all phenotypes have been normalized to have variance of 1.  Under this assumption, a trait's squared [Euclidian norm](https://en.wikipedia.org/wiki/Inner_product_space#Norm_properties) is its heritability: $\lVert v \rVert^2=h^2_v$ where $h^2_v$ is the heritability of $v$.



The above implies that two traits are orthogonal ($\langle u,v \rangle=0$) if and only if they are genetically uncorrelated.



### Perpendicular projection

Let $T_1$ and $T_2$ be the two genetically correlated traits diagrammed as vectors below:

![subtraction-vectors-1](https://github.com/user-attachments/assets/7c084d7b-bb4f-4eac-b118-7d993b5a8e7a)





We aim to decompose $T_1$ into the sum of:

- $F'$, which is perfectly genetically correlated with $T_2$, and
- $R'$, which is orthogonal to (genetically uncorrelated with) $T_2$.

The decomposition is diagrammed below:


![subtraction-vectors-2](https://github.com/user-attachments/assets/a552ee13-08bb-41dd-b7d0-781e4088740b)


Let $P$ denote the perpendicular projector[@halmos1958finite] onto the subspace spanned by $T_2$.  Then

$$
\begin{align}
F'&= PT_1,\\
R'&=(I-P)T_1.
\end{align}
$$

### Interpretation

The primary output of GWAS by subtraction is $R'$, the component of $T_1$ genetically uncorrelated with $T_2$. Studying $R'$ with standard post-GWAS analysis techniques like [MAGMA](MAGMA_Overview.md)[@de2015magma] and [S-LDSC](S_LDSC_For_Cell_And_Tissue_ID.md)[@finucane2018heritability] can shed light on the biological processes important to $T_1$ but absent from $T_2$.


Note that as a linear-algebraic operation, GWAS by subtraction is valid insofar as trait genetics can be approximated by a simple linear model.  While the experience of the last decate and a half of genetics suggest that linear models are very useful, they are necessarily approximations of true biology, which is nonlinear.




## Statistics 


### Joint model

We assume the following data-generating model:

$$
\begin{align}
F &=  x^T \beta_{F} \\
R &=  x^T \beta_{R} \\
T_1 &= \underbrace{a_F F}_{=:F'} + \underbrace{a_R R}_{=:R'} + \delta_1 \label{joint_t_1} \\
T_2&= bF +\delta_2 \label{joint_t_2}  \\
\mathbb{Cov}(F,R)&=0\\
\mathbb{Var}(F)&=1\\
\mathbb{Var}(R)&=1.
\end{align}
$$

Where:


- $T_1,T_2$ are the two traits of interest. We them them as random variables in $\mathbb{R}$.
- There are $M\gg 0$ genetic variants.
- $x\in\mathbb{R}^M$ is the random genotype. We assume $x$ has mean zero, but unlike in [LDSC](LDSC.md), we do not assume it has been variance standardized.  Let $H_i$ be the variance of the $i$th variant.
- $\beta_F,\beta_R\in\mathbb{R}^M$ are the underlying causal effects of the genetic variants.
- $F,R$ are the two orthonormal underlying factors.
- $\delta_1, \delta_2$ are the non-genetic components of the two traits.  We assume these effects are independent of all genotypes.


### Marginal Model


Let's now focus on SNP $i$, and develop a model around the marginal GWAS regression on this SNP.


Define

$$
\begin{align}
\hat\beta_{F,i}&= \frac{\mathrm{Cov}(F, x_i)}{\mathrm{Var}(x_i)}\\
\zeta_{F,i}&=F-\beta_{F,i} x_i\\
\hat\beta_{R,i}&= \frac{\mathrm{Cov}(R, x_i)}{\mathrm{Var}(x_i)}\\
\zeta_{R,i}&=R-\beta_{R,i} x_i\\
\end{align}
$$


This yields the re-written model

$$
\begin{align}
T_1&=a_F F + a_R R + \delta_1\\
T_2& = b F + \delta_2\\
F &= \hat\beta_{F,i}x_i+\zeta_{F,i}\\
R &= \hat\beta_{R,i}x_i+\zeta_{R,i}\\
\mathrm{Cov}(F,R)&=0\\
\mathrm{Cov}(F)&=1\\
\mathrm{Cov}(R)&=1
\end{align}
$$

We assume $\zeta_{F,i},\zeta_{R_i}$ are approximately independent of $x_i$.  This is a good approximation so long as individual variant effects ($\beta_{R,i},\beta_{R,i}$) are small, as is the case for most non-Mendelian traits.

### Theoretical covariance

Next, let us examine the genetic covariance structure of the random variables $(x_i, T_1, T_2)$.  

We will denote by  $\mathrm{GCov}$ and $\mathrm{GVar}$ the genetic covariance and variance respectively\footnote{Because of our earlier assumption that phenotype variance has been normalized to 1, genetic variance equals heritability.}.


$$
\begin{align}
&\mathrm{GCov}(X_i, T_1)\\
&=\mathrm{GCov}(X_i, a_F F + a_R R + \delta_1)\\
&=\mathrm{Cov}(X_i, a_F F + a_R R ) & \text{Since $\delta_1$ is non-genetic}\\
&=\mathrm{Cov}(X_i, a_F (\hat\beta_{F,i}X_i+\zeta_{F,i}) + a_R (\hat\beta_{R,i}X_i+\zeta_{R,i})+\delta_1)\\
&\approx \left(a_F\hat\beta_{F,i}+a_R\hat\beta_{R,i}\right) H_i & \text{By approximate independence}
\end{align}
$$


$$
\begin{align}
&\mathrm{GCov}(X_i, T_2)\\
&=\mathrm{GCov}(X_i, b F+\delta_2)\\
&=\mathrm{Cov}(X_i, b F) & \text{Since $\delta_2$ is non-genetic}\\
&=\mathrm{Cov}(X_i, b (\hat\beta_{F,i}X_i+\zeta_{F,i}))\\
&\approx b\hat\beta_{F,i} H_i  & \text{By approximate independence}
\end{align}
$$


$$
\begin{align}
&\mathrm{GVar}(T_1)\\
&=\mathrm{GVar}(a_F F + a_R R+\delta_1) \\
&=\mathrm{Var}(a_F F + a_R R)& \text{Since $\delta_1$ is non-genetic}\\
&=a_F^2+a_R^2 & \text{Since $F$ and $R$ are uncorrelated}
\end{align}
$$

$$
\begin{align}
&\mathrm{GCov}(T_1,T_2)\\
&= \mathrm{GCov}(a_F F + a_R R+\delta_1, bF+\delta_2)\\
&= \mathrm{Cov}(a_F F + a_R R, bF) & \text{Since $\delta_1,\deta_2$ are non-genetic}\\
&=a_Fb & \text{Since $F$ and $R$ are uncorrelated}
\end{align}
$$


Combining the above yields the following covariance matrix for $(x_i, T_1, T_2)$


$$
\begin{align}
\Sigma_{\text{Theoretical}} &= \begin{bmatrix}
H_i & (a_F\hat\beta_{F_i}+a_R\hat\beta_{R,i})H_i & b\hat\beta_{F,i}H_i \\
(a_F\hat\beta_{F_i}+a_R\hat\beta_{R,i})H_i&  a_F^2+a_R^2 & a_F b   \\
b\hat\beta_{F,i}H_i &a_F b & b^2
\end{bmatrix}
\end{align}
$$

### Empirical covariance

The inputs to GWAS by subtraction are summary statistics for the traits $T_1$ and $T_2$. For SNP i, these will include marginal regression coefficients

$$
\begin{align}
\hat\beta_{T_1,i}&= \frac{\mathrm{Cov}(T_1, x_i)}{H_i}\\
\hat\beta_{T_2,i}&= \frac{\mathrm{Cov}(T_2, x_i)}{H_i}.
\end{align}
$$

Re-arranging, we have





$$
\begin{align}
\mathrm{Cov}(T_1, x_i) &=  \hat\beta_{T_1,i} H_i\\
\mathrm{Cov}(T_2, x_i) &=  \hat\beta_{T_2,i} H_i.
\end{align}
$$

Furthermore, we can apply [LDSC](LDSC.md) and [CT-LDSC](Cross_Trait_LDSC.md) to the $T_1$ and $T_2$ summary statistics to estimate their [genetic covariance](Genetic_Correlation.md) and [heritabilities](Heritability.md) (again, heritability equals genetic variance, since we have assumed that phenotype variances are normalized to 1).  Denote these estimates as $L_{1,2},L_{1,1},L_{2,2}$.


Combining the above, we have that the empirical covariance matrix of $(x_i, T_1, T_2)$ is 

$$
\begin{align}
\Sigma_{\text{Empirical}} &= 
\begin{bmatrix}
H_i & \hat\beta_{T_1,i} H_i & \hat\beta_{T_2,i} H_i\\
\hat\beta_{T_1,i}H_i & L_{1,1} & L_{1,2}\\
\hat\beta_{T_2,i} H_i & L_{1,2} & L_{2,2}
\end{bmatrix}
\end{align}
$$

### Solution


We can equate $\Sigma_{\text{Empirical}}$ and $\Sigma_{\text{Theoretical}}$ to solve for $a_F, a_R, b, \hat\beta_{F,i}, \hat\beta_{R,i}$. We have:


$$
\begin{align}
\Sigma_{\text{Theoretical}} & = \Sigma_{\text{Empirical}}\\
\begin{bmatrix}
H_i & (a_F\hat\beta_{F_i}+a_R\hat\beta_{R,i})H_i & b\hat\beta_{F,i}H_i \\
(a_F\hat\beta_{F_i}+a_R\hat\beta_{R,i})H_i&  a_F^2+a_R^2 & a_F b   \\
b\hat\beta_{F,i}H_i &a_F b & b^2
\end{bmatrix}
&=
\begin{bmatrix}
H_i & \hat\beta_{T_1,i} H_i & \hat\beta_{T_2,i} H_i\\
\hat\beta_{T_1,i}H_i & L_{1,1} & L_{1,2}\\
\hat\beta_{T_2,i} H_i & L_{1,2} & L_{2,2}
\end{bmatrix}
\end{align}
$$

Solving the lower-right $2\times 2$ submatrix, we have:

$$
\begin{align}
b&=\sqrt{L_{2,2}} \label{b_solve} \\
a_F&= \frac{L_{1,2}}{\sqrt{L_{2,2}}} \label{a_F_solve} \\
a_R&=\sqrt{L_{1,1}-\frac{L_{1,2}^2}{L_{2,2}}}  \label{a_R_solve} .
\end{align}
$$


Equating the first columns of the two matrices yields

$$
\begin{align}
\hat\beta_{F,i}&=\frac{\hat\beta_{ T_2,i} }{b} \label{beta_F_solve}\\
\hat\beta_{R,i}&=\frac{1}{a_R}\left(\hat\beta_{ T_1,i} -a_F\frac{\hat\beta_{ T_2,i} }{b}\right) \label{beta_R_solve}.
\end{align}
$$

Note from $(\ref{b_solve}, \ref{a_F_solve}, \ref{a_R_solve})$ that $a_F, a_R$ and $b$ do not depend on the specific genetic variant $i$ under consideration.  This is consistent with the model specified in $(\ref{joint_t_1}, \ref{joint_t_2})$, in which $a_F, a_R$ and $b$ are global constants.


To recap, given summary statistics for traits $T_1$ and $T_2$, we can:

- Run LDSC and CT-LDSC to estimate $L_{1,1},L_{1,2}, L_{2,2}$.
- Apply $(\ref{b_solve},\ref{a_F_solve}, \ref{a_R_solve})$ to estimate $a_F,a_R,$ and $b$.
- Apply $(\ref{beta_F_solve, beta_R_solve})$ to estimate $\hat\beta_{F,i}, \hat\beta_{R,i}$ for each genetic variant $i$.


We would like to pass the $\hat\beta_{R,i}$ to downstream analysis tools like [MAGMA](MAGMA_Overview.md) and [S-LDSC](S_LDSC_For_Cell_And_Tissue_ID.md).  Doing this requires estimates of standard errors.


### Uncertainty

Define $\nu\in\mathbb{R}^5$ to be the key non-redundant entries of $\Sigma_{\text{Empirical}}$.  That is

$$
\begin{align}
\nu:= (\Sigma_{\text{Empirical}, (1,2) },
\Sigma_{\text{Empirical}, (1,3),
\Sigma_{\text{Empirical}, (2,2),
\Sigma_{\text{Empirical}, (2,3),
\Sigma_{\text{Empirical}, (3,3),
)\\
&= (
\hat\beta_{T_1,i}H_i,
\hat\beta_{T_2,i}H_i,
L_{1,1},
L_{1,2},
L_{2,2}
).
\end{align}
$$

Let $\theta\in\mathbb{R}^5$ denote the key parameters we solve for. That is,

$$
\begin{align}
\theta&:= (a_F,a_R,b, \hat\beta_{F,i}, \hat\beta_{R,i})
\end{align}
$$


Let $g:\mathbb{R}^5$



To be continued $\ldots$


