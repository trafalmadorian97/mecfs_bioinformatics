# CT-LDSC


I used [Cross-Trait Linkage Disequilibrium Score Regression](../../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md)[@bulik2015atlas] to compute the genetic correlation between the various LDL GWAS.

Results are below:



{{ include_file("docs/_figs/ldl_ct_ldsc_ct_ldsc_corr_agg_markdown.mdx") }}



The standard error is high.  The genetic correlation of 0.496 is higher than would be expected for two identical traits. This may be in part due to the exclusion of statin users from the Willer GWAS but not the Million Veteran GWAS.