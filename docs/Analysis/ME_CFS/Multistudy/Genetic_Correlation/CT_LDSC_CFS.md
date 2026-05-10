# CT-LDSC


I used [Cross-Trait Linkage Disequilibrium Score Regression](../../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md)[@bulik2015atlas] to compute the genetic correlation between the various CFS-related GWAS, starting with [DecodeME](../../../../Data_Sources/DecodeME.md)[@genetics2025initial] and the [Million Veterans Program](../../../../Data_Sources/Million_Veterans_Program.md)[@verma2024diversity] GWAS of CFS.  Results are below:


{{ include_file("docs/_figs/cfs_ct_ldsc_ct_ldsc_corr_agg_markdown.mdx") }}


The genetic correlation of 0.6523 is high, but not as high as would be expected if the two GWAS were truly measuring an identical trait.  This is consistent with the Million Veterans GWAS using a somewhat imprecise trait definition.  Note also that the standard deviation of the estimated correlation (0.159) is relatively large, perhaps due to the moderate case count in the Million Veterans GWAS.