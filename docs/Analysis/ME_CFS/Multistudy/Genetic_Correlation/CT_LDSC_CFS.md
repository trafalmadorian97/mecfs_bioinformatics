# CT-LDSC


I used [Cross-Trait Linkage Disequilibrium Score Regression](../../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md)[@bulik2015atlas] to compute the [genetic correlation](../../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the various CFS-related GWAS.  I included [DecodeME](../../../../Data_Sources/DecodeME.md)[@genetics2025initial], the subset-GWAS of DecodeME, as well as other GWAS of CFS-like traits.


Results are below:


{{ include_file("docs/_figs/cfs_ct_ldsc_ct_ldsc_corr_agg_markdown.mdx") }}


The first observation is that the various DecodeME subset-GWAS are all highly genetically correlated[^overlap]. Particularly notable is the high genetical correlation between male and female ME/CFS patients; and between ME/CFS patients who do and do not-report an infectious onset. The correlation between male and female patients suggest that the underlying genetic mechanism driving ME/CFS is likely the same in both sexes.  The correlation between infectious-onset and non-infectious onset patients suggests that either: a) the pathological process driving ME/CFS is the same regardless of whether it is triggered by an infectious agent or b) self-reported infectious onset is unreliable.


The genetic correlation between DecodeME and the [Million Veterans Program](../../../../Data_Sources/Million_Veterans_Program.md)[@verma2024diversity] GWAS is 0.6523, which is high, but not as high as would be expected if the two GWAS were truly measuring an identical trait.  This is consistent with the Million Veterans GWAS using a somewhat imprecise trait definition.  Note also that the standard error of the estimated correlation (0.159) is relatively large, perhaps due to the moderate case count in the Million Veterans GWAS.

The genetic correlation between DecodeME and the Neale Lab UK Biobank CFS GWAS is higher, at 0.8017, suggesting a more concordant phenotype.  The increased genetic covariance intercept (0.03582) is consistent with the sample overlap between the two studies induced by DecodeME's use of UK Biobank controls.  Again, the standard error is relatively large, due to the small sample size of the Neale Lab GWAS.


It is interesting that the genetic correlation between the Neale Lab UK Biobank GWAS and the Million Veterans GWAS is lower, at 0.3385.  This may indicate that the two noisy phenotypes are each noisy in a different way. 


Most of the results above suffer from large standard errors.  It may be worthwhile to apply a more precise genetic-correlation estimation method, such as High Definition Likelihood[@ning2020high].


## How to reproduce analysis

Use {{ api_link("this script", "mecfs_bio.analysis.cfs_genetic_corr") }} to reproduce the above analysis.

[^overlap]: This is not a consequence of overlap between some of the subsets.  CT-LDSC automatically accounts for sample overlap using its intercept.