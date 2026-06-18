# GWAS by Subtraction

GWAS by subtraction[@demange2021investigating; @huang2024gwas] is a GenomicSEM[@grotzinger2019genomic] technique that orthogonally decomposes GWAS traits.

Here, we explain GWAS by subtraction twice: once at a high level via linear algebra, and more granularly via statistical modeling.

## Linear algebraic explanation

### Euclidian space

A useful way to understand the GWAS-by-subtraction operation is via linear algebra.


Consider a [Euclidian space](https://en.wikipedia.org/wiki/Euclidean_space) in which:

- GWAS traits are vectors.
- The [inner product](https://en.wikipedia.org/wiki/Inner_product_space) of two traits is their [genetic covariance](Genetic_Correlation.md).  Denote this inner product as $\langle \cdot,\cdot \rangle$.
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


todo