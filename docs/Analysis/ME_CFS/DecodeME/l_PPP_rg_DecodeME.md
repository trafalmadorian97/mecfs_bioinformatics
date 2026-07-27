# PPP CT-LDSC

I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the [DecodeME](../../../Data_Sources/DecodeME.md) GWAS of ME/CFS and the Olink proteomic [@wik2021proximity] GWAS of the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].

## Results

As is standard for LDSC analysis, I restricted the summary statistics to Hapmap 3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. 

[//]: # (T o focus on trans effects, I excluded the cis regions from the proteomic GWAS.)



The results are below:

{{ ppp_rg_data_table("docs/_figs/decode_me_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-decodeme-cis-excluded-rg" )}}



## Interpretation


There are no significant genetic correlations.  Thoughts:

- This could simply be a power issue.  In the future, when we have ME/CFS GWAS with large sample sizes, we may be able to detect genetic correlations with plasma proteins.
- Our [MAGMA](g_MAGMA_DecodeME_Analysis.md) and [S-LDSC](j_S-LDSC_DecodeME_Analysis.md)  analyses suggest that ME/CFS heritability is enriched in neural tissue. Thus, instead of looking for genetic correlations between ME/CFS and plasma protein levels, it could be more fruitful to look for genetic correlations between ME/CFS and cerebro-spinal fluid protein levels.



