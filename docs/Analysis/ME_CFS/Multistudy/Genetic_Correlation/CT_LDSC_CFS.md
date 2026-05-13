# CT-LDSC


I used [Cross-Trait Linkage Disequilibrium Score Regression](../../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md)[@bulik2015atlas] to compute the genetic correlation between the various CFS-related GWAS.


Results are below:


{{ include_file("docs/_figs/cfs_ct_ldsc_ct_ldsc_corr_agg_markdown.mdx") }}


The genetic correlation between [DecodeME](../../../../Data_Sources/DecodeME.md)[@genetics2025initial] and the [Million Veterans Program](../../../../Data_Sources/Million_Veterans_Program.md)[@verma2024diversity] GWAS is 0.6523, which is high, but not as high as would be expected if the two GWAS were truly measuring an identical trait.  This is consistent with the Million Veterans GWAS using a somewhat imprecise trait definition.  Note also that the standard error of the estimated correlation (0.159) is relatively large, perhaps due to the moderate case count in the Million Veterans GWAS.


The genetic correlation between DecodeME and the Neale Lab UK Biobank CFS GWAS is higher, at 0.8017, suggesting a more concordant phenotype.  The increased genetic covariance intercept (0.03582) is consistent with the sample overlap between the two studies induced by DecodeME's use of UK Biobank controls.  Again, the standard error is relatively large, due to the small sample size of the Neale Lab GWAS.


It is interesting that the genetic correlation between the Neale Lab UK Biobank GWAS and the Million Veterans GWAS is lower, at 0.3385.  This may indicate that the two noisy phenotypes are each noisy in a different way. 


All of the results above suffer from large standard errors.  It may be worthwhile to apply a more precise genetic-correlation estimation method, such as High Definition Likelihood[@ning2020high].

