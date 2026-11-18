---
tags:
  - ABF colocalization
---
# PPP Coloc
For further insights into the molecular mechanisms of schizophrenia, I ran [ABF Colocalization](../../../Bioinformatics_Concepts/Colocalization.md)[@giambartolomei2014bayesian] to search for loci at which the 2022 PGC GWAS of Schizophrenia[@trubetskoy2022mapping] colocalized with proteomics GWAS from the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md)[@sun2023plasma].

All locus/protein pairs with greater than 80 percent probability of colocalization are shown below


{{
coloc_results_table("docs/_figs/pgc_2022_sch_ppp_coloc_full_ppp_coloc_display.parquet","pgc_sch_ppp_coloc")
}}

There are a large number of protein/locus pairs with a high probability of colocalization.

Some interesting observations:

- There are three colocalizations with cis pQTLs: CD40, RABEP1, and PLEKHO1. [CD40](https://www.genecards.org/card/CD40) has a known role in the immune system, [RABEP1](https://www.genecards.org/card/RABEP1) is involved in vesicular traffic, while [PLEKHO1](https://www.genecards.org/card/PLEKHO1) is involved in actin dynamics

- ICAM5, the only protein found to Bonferroni-significantly genetically correlated with schizophrenia in my [CT-LDSC analysis](SC_PGC_PPP_RG.md), did not colocalize.  This suggests a global polygenic overlap of genetic architecture, but no sufficiently strong local overlap at any particular locus.