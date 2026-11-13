---
tags:
  - ABF colocalization
---
# Exploratory PPP Coloc

I ran [ABF Colocalization](../../../Bioinformatics_Concepts/Colocalization.md)[@giambartolomei2014bayesian] to search for loci at which [DecodeME](../../../Data_Sources/DecodeME.md) GWAS 1[@genetics2025initial] colocalized with Olink[@wik2021proximity] proteomic GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].

This was an exploratory analysis, so I included all DecodeME loci at which the lead DecodeME variant had a p value less that $10%^{-5}$.  Note that this is a looser threshold than what is normally used for genome-wide significance.  The main results are below[^trim_note]:


{{ coloc_results_table("docs/_figs/decode_me_gwas_1_exploratory_2_ppp_coloc_display.parquet", id="ukbb-ppp-coloc-decodeme" )}}


None of the genome-wide significant DecodeME loci colocalized with a proteomic GWAS.  Two sub-genome-wide-significant DecodeME loci did colocalize.


Below, I show the DecodeME Manhattan plot at these two loci, together the Manhattan plots of the colocalizing proteomic GWAS:


{{
static_img_embed("docs/_figs/decode_me_gwas_1_exploratory_2_ppp_coloc_stack_plots_png/chr1_7941726-8941726.png")
}}


{{
static_img_embed("docs/_figs/decode_me_gwas_1_exploratory_2_ppp_coloc_stack_plots_png/chr12_123424955-124424955.png")
}}



I looked up the colocalizing proteins in some proteomics databases.  The results follow.


{{
data_table("docs/_figs/decode_me_gwas_1_exploratory_2_ppp_coloc_protein_info.parquet",id="decode_me_ppp_coloc_info")
}}




[^trim_note]: The table is trimmed to include only proteins with $P(H_4)\ge 0.8$.