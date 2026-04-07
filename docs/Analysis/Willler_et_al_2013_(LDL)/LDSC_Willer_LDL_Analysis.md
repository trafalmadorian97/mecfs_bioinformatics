# LDSC

I applied [Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to the [LDL](../Verma_et_al_(LDL)/LDL_Overview.md) GWAS of Willer et al.[@global2013discovery].  The aim was to estimate heritability and look for signs of stratification and confounding.


The results are below:



--8<--  "docs/_figs/willer_et_al_2023_ldl_eur_ldsc_heritability_markdown.mdx"


In is interesting to contrast these results with those of the [Million Veterans Program GWAS of LDL](../Verma_et_al_(LDL)/LDSC_LDL_Analysis.md).  The most notable difference is in the LDSC intercept.  While the Million Veterans GWAS yields an LDSC intercept of 1.586, which strongly suggests population stratification or confounding, the Willer et al. GWAS yields an intercept of 0.9984, which does not suggest stratification or confounding. It is unclear exactly what explains this difference, but it does indicate that the Willer GWAS should be treated as more trustworthy than the Million Veterans GWAS.



The other key difference is in heritability.  The heritability of 0.1922 found by the Willer et al. GWAS is significantly higher than the heritability of 0.1145 found by the Million Veterans GWAS.  Given that the Million veterans intercept indicates stratification or confounding, I put more trust in the heritability estimate from Willer et al.