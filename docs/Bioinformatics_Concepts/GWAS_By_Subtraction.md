# GWAS By Subtraction

GWAS by subtraction[@demange2021investigating; @huang2024gwas] is a technique based on GenomicSEM[@grotzinger2019genomic] that orthogonally decomposes GWAS traits.

Here, we explain GWAS by subtraction twice: once at a high-level via linear algebra, and once a more granular level of detail via explicit probabilistic modeling.

## Linear algebraic explanation

### The Euclidian space of GWAS traits

A useful way to understand the GWAS-by-subtraction operation is via linear algebra.


Consider a [Euclidian space](https://en.wikipedia.org/wiki/Euclidean_space) in which:

- GWAS traits are vectors.
- The [inner product](https://en.wikipedia.org/wiki/Inner_product_space) of two traits is their [genetic covariance](Genetic_Correlation.md).  Denote this inner product as $\langle \cdot,\cdot \rangle$.
- We assume all phenotypes have been normalized to have variance of 1.  Under this assumption, the squared [Euclidian norm](https://en.wikipedia.org/wiki/Inner_product_space#Norm_properties) of a trait is that trait's heritability: $\lVert v \rVert^2=h^2_v$ where $h^2_v$ is the heritability of $v$.



From the above, it further follows that two traits are orthogonal ($\langle u,v \rangle=0$) if and only if they are genetically uncorrelated.



### GWAS by subtraction as perpendicular projection

Now suppose we have two GWAS traits $T_1, T_2$.  We aim to decompose $T_1$ into:

- A component $F'$ that is perfectly genetic correlated with $T_2$. 
- A component $R'$ that is orthogonal (genetically uncorrelated) with $T_2$.





## Statistical explanation