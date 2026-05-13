# LDSC Analysis

I applied [Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to the non-23-and-me European subset of the adult height GWAS of Yengo et al.[@yengo2022saturated] to estimate [heritability](../../../Bioinformatics_Concepts/Heritability.md) and look for signs of population stratification or confounding.

The results are below:



{{ include_file("docs/_figs/yengo_et_al_height_ldsc_heritability_markdown.mdx") }}


A SNP heritability of approximately 43% is consistent with estimates from other sources.

An LD-score regression intercept of 2.429 and a genomic control $\lambda$ of 5.191 may seem high.  Note however, that height is a highly polygenic trait, and the Yengo GWAS has a very large sample size.  This is reflected in the very large mean $\chi^2$ (15.38).  In such a circumstances, the best measure of degree of potential stratification/ confounding is the attenuation ratio, which is calculated as 

$$
\frac{\text{intercept} -1}{\chi^2_{\text{mean}} -1}.
$$

Here the attenuation ratio is 0.09943, which suggests an acceptably low level of stratification.