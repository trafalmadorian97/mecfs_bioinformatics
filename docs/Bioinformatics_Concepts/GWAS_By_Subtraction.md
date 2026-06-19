# GWAS by Subtraction

GWAS by subtraction[@demange2021investigating; @huang2024gwas] is a GenomicSEM[@grotzinger2019genomic] technique that orthogonally decomposes GWAS traits.

Here, we explain GWAS by subtraction twice: once at a high level via linear algebra, and more granularly via statistical modeling.

## Linear algebraic explanation

### Euclidian space

A useful way to understand the GWAS-by-subtraction operation is via linear algebra.


Consider a [Euclidian space](https://en.wikipedia.org/wiki/Euclidean_space) in which:

- GWAS traits are vectors.
- The [inner product](https://en.wikipedia.org/wiki/Inner_product_space) of two traits is their [genetic covariance](Genetic_Correlation.md).  Denote the inner product of $u$ and $v$ as $\langle u,v \rangle$.
- We assume all phenotypes have been normalized to have variance of 1.  Under this assumption, the squared [Euclidian norm](https://en.wikipedia.org/wiki/Inner_product_space#Norm_properties) of a trait is that trait's heritability: $\lVert v \rVert^2=h^2_v$ where $h^2_v$ is the heritability of $v$.



From the above, it further follows that two traits are orthogonal ($\langle u,v \rangle=0$) if and only if they are genetically uncorrelated.



### Perpendicular projection

Suppose we have two genetically correlated GWAS traits $T_1, T_2$. These traits are diagrammed as vectors below

![subtraction-vectors-1](https://github.com/user-attachments/assets/7c084d7b-bb4f-4eac-b118-7d993b5a8e7a)





We aim to decompose $T_1$ into the sum of:

- $F'$, which is perfectly genetically correlated with $T_2$, and
- $R'$, which is orthogonal (genetically uncorrelated) to $T_2$.

The orthogonal vector sum is diagrammed below:


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


Note that as a linear-algebraic operation, GWAS by subtraction is valid insofar as trait genetics can be approximated by a simple linear model.  While the experience of the last decate and a half of genetics suggest that linear models are very useful, they are necessarily approximations of true biology, which is inherently nonlinear.




## Statistical explanation


### Joint model

We assume the following data-generating model

$$
\begin{align}
F &=  x^T \beta_{F} + \epsilon_F\\
R &=  x^T \beta_{R} + \epsilon_R\\
T_1 &= \underbrace{a_F F}_{=:F'} + \underbrace{a_R R}_{=:R'}\\
T_2&= bF\\
\mathbb{Cov}(F,R)&=0\\
\mathbb{Var}(F)&=1\\
\mathbb{Var}(R)&=1
\end{align}
$$

Where:
- $T_1,T_2$ are the two traits of interest. We them them as random variables in $\mathbb{R}$.
- There are $M\gg 0$ genetic variants.
- $x\in\mathbb{R}^M$ is the random genotype. We assume $x$ has mean zero, but unlike in [LDSC](LDSC.md), we do not assume it has been variance standardized.  
- $\beta_F,\beta_R\in\mathbb{R}^M$ are the underlying causal effects of the genetic variants.
- $F,R$ are the two orthonormal underlying factors.


### Marginal Model


Let's now focus on SNP $i$, and develop a model around the marginal GWAS regression on this SNP.


Define

$$
\begin{align}
\hat\beta_{F,i}&= \frac{\mathrm{Cov}(F, x_i)}{\mathrm{Var}(x_i)}\\
\zeta_{F,i}&=F-\hat\beta_{F,i} x_i\\
\hat\beta_{R,i}&= \frac{\mathrm{Cov}(R, x_i)}{\mathrm{Var}(x_i)}\\
\zeta_{R,i}&=R-\hat\beta_{R,i} x_i\\
\end{align}
$$


This yields the re-written model

$$
\begin{align}
T_1&=a_F F + a_R R\\
T_2& = b F\\
F&= \hat\beta_{F,i}x_i+\zeta_{F,i}\\
R &= \hat\beta_{R,i}x_i+\zeta_{R,i}\\
\mathrm{Cov}(F,R)&=0\\
\mathrm{Cov}(F)&=1\\
\mathrm{Cov}(R)&=1
\end{align}
$$

### Theoretical covariance matrix

Next, let us examine the covariance structure of the random variables $(x_i, T_1, T_2)$.  


To be continued $\ldots$

