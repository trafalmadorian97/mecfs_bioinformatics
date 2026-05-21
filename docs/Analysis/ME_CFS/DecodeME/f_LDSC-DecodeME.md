# LDSC

I applied [Linkage Disequlibirum Score Regression](../../../Bioinformatics_Concepts/LDSC.md)[@bulik2015ld] to GWAS-1 from [DecodeME](../../../Data_Sources/DecodeME.md) to estimate heritability and look for evidence of stratification.

The results are in the table below:


{{include_file("docs/_figs/decode_me_gwas_1_ldsc_heritability_markdown.mdx")}}


A liability scale heritability of 0.0814 is similar but not identical to the value reported in the original DecodeME preprint[@genetics2025initial].  This discrepancy may perhaps be explained by slight differences in the technique used to liftover from genome build 38 to 37.

Reassuringly, the attenuation ratio is less than 1, indicating no evidence of confounding due to population stratification.