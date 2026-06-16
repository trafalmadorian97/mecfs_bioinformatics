# Overview

I used GWAS-by-subtraction[@demange2021investigating; @huang2024gwas] to decompose the ME/CFS trait measured by [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] into the sum of two traits:

- A factor trait $F$ that is perfectly [genetically correlated](../../../../../Bioinformatics_Concepts/Genetic_Correlation.md) with the trait from Johnston et al.'s GWAS of multisite pain[@johnston2019genome].
- A remainder trait $R$, that is genetically uncorrelated with the trait from Johnston et al.'s GWAS.


I then analyzed the remainder trait $R$ using my standard suite for post-GWAS analyses.   Previously, we have observed both a [high genetic correlation](../../../Genetic_Correlation/CT-LDSC_Study.md) between multi-site pain and ME/CFS, and a high degree of concordance between local significance patterns of the two GWAS.  The aim of the present investigation is to separately analyze $R$, that part of the DecodeME GWAS signal which is distinct from the multisite pain GWAS signal.  Clearly, the symptoms of multisite pain and ME/CFS are different. I hoped to find the genetics underlying this difference.

Mathematically, GWAS-by-subtraction is a purely linear-algebraic operation: we define a Euclidian space in which GWAS traits are vectors and inner product between two vectors is their [genetic covariance](../../../../../Bioinformatics_Concepts/Genetic_Correlation.md).  In this Euclidian space, $R$ is the perpendicular projection of ME/CFS onto the orthogonal complement of chronic pain[^algebra_note].


[^algebra_note]:  For more on the underlying mathematics, see Halmos's classic text[@halmos1958finite].