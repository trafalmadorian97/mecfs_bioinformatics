---
tags:
  - CT-LDSC
---
# PPP CT-LDSC

I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the Bellenguez et al.[@bellenguez2022new] meta-GWAS of Alzheimer's disease and Olink proteomic [@wik2021proximity] GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].

## Results

As is standard for LDSC analysis, I restricted the summary statistics to Hapmap[@international2005haplotype] 3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. To focus on trans effects, I excluded the cis regions from the proteomic GWAS.



The results are below:

{{ ppp_rg_data_table("docs/_figs/bellenguez_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-alz-rg" )}}


## Interpretation

 Unlike some of the [other PPP genetic correlation analyses](../../RA/DECODE_(Seropositive)/6_RA_Seropositive_ppp_rg.md), I was not able to find strong biological or medical literature supporting the importance of the top most significant proteins.  WFIKKN2, the top hit, has some minor supporting evidence: a Singaporean proteomic study of memory-clinic patients found in a secondary analysis that WFIKKN2 was associated with changes in cognitive function over time across all cognitive subgroups[@sim2025plasma].
