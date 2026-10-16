# LDSC


[Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] was applied to the Pulit et al.[@pulit2019meta] meta-analysis of BMI.


The aim was to estimate heritability and look for signs of stratification and confounding.

The results are below:


{{markdown_table("docs/_figs/pulit_et_al_2018_bmi_ldsc_heritability_markdown.mdx")}}


To illustrate the fit of the LDSC model the BMI GWAS, I generated the following LDSC diagnostic plot:


{{
ldsc_diagnostic_plot(src="docs/_figs/pulit_et_al_2018_bmi_ldsc_diagnostic_plot/ldsc_diagnostic.html", id="bmi_ldsc")
}}


LDSC produces a heritability estimate of approximately $0.2$. While the LDSC intercept is significantly greater than 1, the [attenuation ratio ](../../../Bioinformatics_Concepts/LDSC.md#the-attenuation-ratio) is small, indicating relatively low risk of confounding.
