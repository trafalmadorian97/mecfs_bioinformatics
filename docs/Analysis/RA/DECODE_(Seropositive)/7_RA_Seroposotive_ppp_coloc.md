---
tags:
- ABF colocalization
---
# PPP Coloc
For further insights into the mechanisms of seropositive rheumatoid arthritis, I ran [ABF Colocalization](../../../Bioinformatics_Concepts/Colocalization.md)[@giambartolomei2014bayesian] to search for loci at which the DECODE meta-GWAS of seropositive RA[@saevarsdottir2022multiomics] colocalized with proteomics GWAS from the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md)[@sun2023plasma]. 


The results are below.


{{
coloc_results_table("docs/_figs/decode_ra_seropositive_2000_ppp_coloc_display.parquet","ra_ppp_coloc")
}}



There are a large number of protein/locus pairs with a high probability of colocalization.

Some interesting observations

- Consistent with the known pathobiology of rheumatoid arthritis, many of these proteins have known immunological roles.
- There are three proteins that colocalize with RA at cis-pQTL: CD40, FCGR2A, SPAG1. Of these, CD40 and FCGR2A have clear immune functions.
- Looking up CCL19, the protein I found to be most significantly genetically correlated with seropositive RA in my  [CT-LDSC study](6_RA_Seropositive_ppp_rg.md), we see that the CCL19 proteomic GWAS colocalizes with RA at three separate loci _chr1:112716068-114261186_, _chr6:158569915-15956991_, and _chr2:203380280-204380280_.  This reinforces the observation that there is a high degree of genetic concordance between the genetic architectures of CCL19 and rheumatoid arthritis.