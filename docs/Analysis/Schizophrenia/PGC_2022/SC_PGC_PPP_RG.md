---
tags:
  - CT-LDSC
---
# PPP CT-LDSC

I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the 2022 PGC GWAS of schizophrenia[@trubetskoy2022mapping] and the Olink proteomic [@wik2021proximity] GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].

## Results

As is standard for LDSC analysis, I restricted the summary statistics to Hapmap[@international2005haplotype] 3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. To focus on trans effects, I excluded the cis regions from the proteomic GWAS.


The results are below:



{{ ppp_rg_data_table("docs/_figs/pgc2022_sch_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-pgc-sch-rg" )}}

## Interpretation

ICAM5 is the only Bonferroni-significant protein.  [The NCBI description of ICAM5](https://www.genecards.org/card/ICAM5) says that _"This protein is expressed on the surface of telencephalic neurons and displays two types of adhesion activity, homophilic binding between neurons and heterophilic binding between neurons and leukocytes. It may be a critical component in neuron-microglial cell interactions in the course of normal development or as part of neurodegenerative diseases"_.  Thus, at a high level, it is at least plausible that ICAM5 should play a role in schizophrenia, a disease of the central nervous system.  It is unclear, however, how increased levels of ICAM5 in the plasma would relate to pathological processes in the central nervous system, given that in theory the blood-brain barrier should separate the CNS from the plasma.
