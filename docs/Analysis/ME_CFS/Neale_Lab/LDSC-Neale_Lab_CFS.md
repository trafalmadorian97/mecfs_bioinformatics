# LDSC

I applied [Linkage Disequlibirum Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to the Neale Lab's GWAS of self-reported chronic fatigue syndrome (CFS) in the UK Biobank.

Since the Neale Lab's CFS phenotype is derived from a single-question self-report, it is less precise than DecodeME's detailed phenotyping.  This imprecision combined with the relatively small number of cases is expected to produce noisy summary statistics.  Nevertheless, analyzing these summary statistics is still interesting and worthwhile.


The LDSC results are shown below:




{{ include_file("docs/_figs/neale_lab_cfs_ldsc_heritability_markdown.mdx") }}

The liability-scale heritability of 8% is similar to that [seen in DecodeME](../DecodeME/f_LDSC-DecodeME.md).

The low intercept and moderate attenuation ratio suggest low risk of confounding by stratification.
