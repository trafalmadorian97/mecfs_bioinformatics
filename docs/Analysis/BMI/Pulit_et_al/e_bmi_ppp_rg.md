---
tags:
  - CT-LDSC
---
# PPP CT-LDSC
I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the Pullit et al.[@pulit2019meta] meta-GWAS of body mass index (BMI) and proteomic GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].


Note that because of its assumptions of uniform polygenicity, CT-LDSC measures genetic correlation due to diffuse polygenic effects.  It does not measure genetic correlation due to highly-concentrated locus-specific effects.


## Results

As is standard for LDSC analysis, I restricted to the statistics to Hapmap3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. Because my [previous heritability experiment](../../Proteomics/UKBB_PPP/PPP_LDSC.md) suggested little difference between using all SNPs and excluding the cis-region near the protein of interest, I only ran this experiment with the cis region excluded.

The results are below:

{{ ppp_rg_data_table("docs/_figs/pulit_et_al_2018_bmi_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-bmi-rg" )}}

## Interpretation

There are large number of significantly genetically correlated proteins, many of which are consistent with the existing literature

The most strongly genetically correlated protein is leptin (LEP).  This protein is known to produced by fat cells in proportion to fat mass[@yeo2018gene], so its presence make sense


