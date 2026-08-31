---
tags:
  - LDSC
---
# LDSC

I applied [Linkage Disequlibirum Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to [DecodeME](../../../Data_Sources/DecodeME.md) GWAS-1 to estimate heritability and look for evidence of stratification.

The results follow:


{{include_file("docs/_figs/decode_me_gwas_1_ldsc_heritability_markdown.mdx")}}


A liability scale heritability of 0.0814 is similar but not identical to the value reported in the original DecodeME preprint[@genetics2025initial].  This discrepancy may be explained by slight differences in the technique used to liftover from genome build 38 to 37.

The LDSC intercept is less than one. On the one hand, this is reassuring in the sense that it provides no evidence of confounding due to population stratification. On the other hand, an intercept less than 1 is a mild violation of assumptions underlying the LDSC model, and thus is a diagnostic flag that should be followed-up on.