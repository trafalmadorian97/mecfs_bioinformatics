# PPP CT-LDSC
I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the [Million Veterans Program](../../../Data_Sources/Million_Veterans_Program.md)[@verma2024diversity] GWAS of myocardial infarction and Olink proteomic assay[@wik2021proximity] GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].


## Results

As is standard for LDSC analysis, I restricted the summary statistics to Hapmap 3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. Because my [previous heritability experiment](../../Proteomics/UKBB_PPP/PPP_LDSC.md) suggested little difference between using all SNPs and excluding the cis-region near the protein of interest, I only ran this experiment with the cis region excluded.

The results are below:

{{ ppp_rg_data_table("docs/_figs/mv_myocardial_infarction_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-mv-mi-rg" )}}



## Interpretation


In contrast to some of my previous genetic correlation analysis against UKBB PPP data, here we see a very large number of Bonferroni-significant proteins.

The protein with the most significant genetic correlation with myocardial infarction is GDF-15 (Growth differentiation factor 15).  This finding is consistent with known biology.  For instance, Kato et al.[@kato2023growth] found in a meta-analysis of data from 8 clinical trials that GDF-15 levels were strongly predictive of future myocardial infarction and other adverse cardiac events among patients with stable atherosclerotic cardiovascular disease or patients stabilized after acute coronary syndrome.  See the graphical abstract from this paper below:

![kato-abstract](https://github.com/user-attachments/assets/29b267f9-c9a4-4db6-a965-46995027cd36)



## Caveats

Recall that while CT-LDSC is generally not vulnerable to environmental confounding, it tells us nothing about the causal direction between myocardial infarction and proteins of interest.


## How to reproduce

To reproduce these results, run the script  {{api_link("here","mecfs_bio.analysis.myocardial_infarction_ppp_rg_analysis")}}.