 # LDSC

I applied [Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to the [Million Veterans Program](../../Data_Sources/Million_Veterans_Program.md)[@verma2024diversity] migraine GWAS.  The aim was to estimate heritability and look for signs of stratification and confounding.

The results are below:



{{ include_file("docs/_figs/million_veterans_migraine_eur_ldsc_heritability_markdown.mdx.") }}


A heritability of 0.1005 is low-to-moderate.

However, an attenuation ratio of 0.3561 is quite high, suggesting possible confounding by population stratification.  This suggests we should be wary of results generated from analysis of these summary statistics.  While LDSC-based methods are likely trustworthy, since they model and remove the stratification component, non-LDSC methods may be less trustworthy and could produce misleading results.





