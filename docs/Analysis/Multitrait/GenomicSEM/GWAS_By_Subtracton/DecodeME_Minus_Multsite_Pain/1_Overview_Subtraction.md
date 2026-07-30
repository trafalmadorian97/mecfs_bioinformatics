# Overview

I used [GWAS by subtraction](../../../../../Bioinformatics_Concepts/GWAS_By_Subtraction.md)[@demange2021investigating; @huang2024gwas] to decompose the [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] GWAS of ME/CFS into the sum of two traits[^algebra_note]:

- A factor trait $F$ that is perfectly [genetically correlated](../../../../../Bioinformatics_Concepts/Genetic_Correlation.md) with Johnston et al.'s GWAS of multisite pain[@johnston2019genome].
- A remainder trait $R$, that is genetically uncorrelated with the trait from Johnston et al.'s GWAS.


I then analyzed the remainder trait $R$ using my standard suite for post-GWAS analyses.   Previously, [we observed](../../../Genetic_Correlation/CT-LDSC_Study.md) a high [genetic correlation](../../../../../Bioinformatics_Concepts/Genetic_Correlation.md) between multisite pain and ME/CFS, and a high degree of concordance between local significance patterns of the two GWAS.  Despite this genetic correlation, it is clear that ME/CFS and multsite pain are distinct conditions. By studying $R$, that part of ME/CFS that is independent of multisite pain, I aimed to find genetics underlying the distinction between the two conditions.



[^algebra_note]: Mathematically, GWAS-by-subtraction is a purely linear-algebraic operation: we define a Euclidian space in which GWAS traits are vectors and inner product between two vectors is their [genetic covariance](../../../../../Bioinformatics_Concepts/Genetic_Correlation.md).  In this Euclidian space, $R$ is the perpendicular projection of ME/CFS onto the orthogonal complement of chronic pain. For more on the underlying mathematics, see Section 75 _"Perpendicular Projection"_ in Halmos's classic linear algebra text[@halmos1958finite].