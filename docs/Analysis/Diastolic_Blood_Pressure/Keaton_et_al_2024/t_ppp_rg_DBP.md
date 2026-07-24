# PPP Genetic Correlation

I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between Keaton et al.'s meta-GWAS of diastolic blood pressure (DBP)[@keaton2024genome] and GWAS of Olink proteomic assays[@wik2021proximity] from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].

## Results

As is standard for LDSC analysis, I restricted to the statistics to Hapmap3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. Because my [previous heritability experiment](../../Proteomics/UKBB_PPP/PPP_LDSC.md) suggested little difference between using all SNPs and excluding the cis-region near the protein of interest, I only ran this experiment with the cis region excluded.

The results are below:

{{ ppp_rg_data_table("docs/_figs/keaton_dbp_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-dbp-rg" )}}


## Interpretation


There is a single protein whose genetic correlation with DBP is Bonferroni-significant, and multiple proteins with genetic correlations significant according to the Benjamini-Hochberg[@benjamini1995controlling] procedure.


CTSB (Cathepsin B), the Bonferroni-significant protein, makes biological sense: in a mouse model of [nephrotic syndrome](https://en.wikipedia.org/wiki/Nephrotic_syndrome), inhibition of CTSB prevented hypertension[@larionov2019cathepsin].