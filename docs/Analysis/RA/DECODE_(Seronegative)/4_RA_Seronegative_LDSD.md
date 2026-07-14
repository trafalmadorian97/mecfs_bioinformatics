# LDSC

[Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] was applied to the DECODE meta-GWAS of seronegative rheumatoid arthritis[@saevarsdottir2022multiomics].  The aim was to estimate heritability and look for signs of stratification and confounding.

The results are below:


{{markdown_table("docs/_figs/decode_ra_seronegative_ldsc_heritability_markdown.mdx")}}

The liability scale heritability is around 0.05, approximately half the [liability scale heritability of seropositive RA](../DECODE_(Seropositive)/4_RA_Seropositive_LDSC.md).  Note however, that this finding is dependent on our estimates of the population prevalence of the two types of RA.

The relatively low intercept suggests a relatively low risk of confounding by population stratification.