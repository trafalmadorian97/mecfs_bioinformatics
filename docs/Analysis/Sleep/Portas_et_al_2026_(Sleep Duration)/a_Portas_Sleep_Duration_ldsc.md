---
tags:
- LDSC
---
# LDSC

I applied [Linkage Disequlibirum Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to Portas et al.'s GWAS of sleep duration[@portas2026genetic] to estimate heritability and look for evidence of [stratification](../../../Bioinformatics_Concepts/Population_Stratification.md).



The results follow:


{{include_file("docs/_figs/portas_et_al_2026_sleep_duration_ldsc_heritability_markdown.mdx")}}


The intercept and attention ratio suggest a low degree of confounding due to population stratification. A heritability of about 12 percent is low to moderate.


The diagnostic plot below indicates an overall good match between the GWAS data and the LDSC model of polygenic heritability.



{{
ldsc_diagnostic_plot(src="docs/_figs/portas_et_al_2026_sleep_duration_ldsc_diagnostic_plot/ldsc_diagnostic.html",id="sleep_duration_ldsc_diagnostic")
}}




