# LDSC 

I applied [Linkage Disequlibirum Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to the GWAS of chronic fatigue syndrome (CFS) performed as part of the [Million Veterans Program](../../../Data_Sources/Million_Veterans_Program.md)[@verma2024diversity].  Note that while [DecodeME](../../../Data_Sources/DecodeME.md)[@genetics2025initial] used a carefully constructed case definition, Million Veterans GWAS are based on PheCodes, which are aggregations of billing codes.  We can thus expect that the Million Veterans GWAS may be significantly more noisy than DecodeME.

The results of running LDSC are listed below:


{{ include_file("docs/_figs/million_veterans_cfs_ldsc_heritability_markdown.mdx") }}

Notably, the liability scale heritability is less than half that found in the analogous DecodeME analysis.

While the LDSC intercept is low, so is the mean $\chi^2$ statistic.  Consequently, the attenuation ratio is moderate-to-high.  This suggests some possibility of confounding by stratification.
