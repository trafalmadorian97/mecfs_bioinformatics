# LD Score Reference Data
An LD score measures the extent to which a SNP is correlated with other neighbouring SNPs in the context of a given population.
The concept of the LD score was introduced [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/) by Bulik-Sullivan et al, and has since been widely used in the computation of heritability, genetic correlation, and other important quantities.  

Formally, $l_i$, the LD score for SNP $i$ defined as

$$l_i := \sum_j r_{i,j}^2 $$

where $r_{i,j}$ is the Pearson correlation coefficient between SNP $i$ and SNP $j$ in the population of interest.