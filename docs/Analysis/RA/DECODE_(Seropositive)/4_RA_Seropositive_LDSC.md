# LDSC

[Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] was applied to the DECODE meta-GWAS of seropositive rheumatoid arthritis[@saevarsdottir2022multiomics].  The aim was to estimate heritability and look for signs of stratification and confounding.

The results are below:


{{markdown_table("docs/_figs/decode_ra_seropositive_ldsc_heritability_markdown.mdx")}}

We get a liability-scale heritability of 0.1003, and an intercept indicating a relatively low risk of population stratification.

The original DECODE paper reports an observed-scale heritability of 0.19. To check whether our heritability result is consistent with this, we can convert from liability- to observed-scale heritability.  Let $h^2_O$ denote observed-scale heritability, $h^2_L$  liability-scale heritability, $K$ population prevalence, $P$ sample prevalence, $\Phi^{-1}$ the inverse of the normal cumulative distribution function, and $\phi$ the normal density function.

We apply the [standard formula](https://cloufield.github.io/gwaslab/HeritabilityConversion/) and substitute the same prevalence parameters we used {{api_link(text="in our LDSC run",module_path="mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_standard_analysis")}}.



$$
\begin{align}
h^2_O&= h^2_L \frac{P(1-P)  }{ K^2 (1-K)^2   } \phi(\Phi^{-1}(1-K))^2\\
&=0.1003 \frac{0.5^2}{0.005(1-0.005)} \phi(\Phi^{-1}(1-0.005))^2\\
&\approx .1003 (9900) (0.000209)\\
&\approx .1003 (2.07)\\
&\approx 0.21
\end{align}
$$

As expected, this result approximates the observed-scale heritability from the DECODE paper (0.19).  The remaining discrepancy is likely due to differences in SNP filtering.