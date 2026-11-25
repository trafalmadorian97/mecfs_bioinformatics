---
tags:
  - LDSC
---

# LDSC
I applied [Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to the  CRP GWAS of Said et al.[@verma2024diversity] to estimate [heritability](../../../Bioinformatics_Concepts/Heritability.md) and look for signs of population stratification or confounding.


The results are in the table below:


{{ include_file("docs/_figs/said_et_al_2022_crp_eur_ldsc_heritability_markdown.mdx") }}


A heritability of 0.1683 is low-to-moderate, while a LDSC intercept 1.065 suggested a well-structured GWAS without obvious confounding or stratification. The high mean chi-squared (2.813) indicates a highly polygenic trait.

To illustrate the fit of the LDSC model to the data, I also generated the LDSC diagnostic plot below.


{{
ldsc_diagnostic_plot(src="docs/_figs/said_et_al_2022_crp_eur_ldsc_diagnostic_plot/ldsc_diagnostic.html", id="said_ldsc")
}}

