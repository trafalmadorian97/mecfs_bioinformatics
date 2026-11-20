---
tags:
  - ABF colocalization
---
# PPP Coloc

For further insight into the molecular mechanisms of Alzheimer's I ran [ABF Colocalization](../../../Bioinformatics_Concepts/Colocalization.md)[@giambartolomei2014bayesian] to search for loci at which the Bellenguez et al.[@bellenguez2022new] GWAS of Alzheimer's colocalized with proteomics GWAS from the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md)[@sun2023plasma].



All locus/protein pairs with greater than 80 percent probability of colocalization are listed below


{{
coloc_results_table("docs/_figs/bellenguez_alzheimers_ppp_coloc_full_ppp_coloc_display.parquet","alz_ppp_coloc")
}}


- There are a significant number of high probability colocalizations.
- Notably, there is a strong cis-colocalization with GRN (granulin precursor).  This is in agreement with the [known association of certain GRN variants with frontotemporal dementia.](https://www.ncbi.nlm.nih.gov/gene/2896).
- It is also notable there are number of colocalizations with complement-related genes such as CR1, C1S, and C1R.



