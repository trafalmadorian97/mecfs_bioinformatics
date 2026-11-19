# PPP Coloc
For further insights into the molecular mechanisms of Fibromyalgia, I ran [ABF Colocalization](../../../Bioinformatics_Concepts/Colocalization.md)[@giambartolomei2014bayesian] to search for loci at which the Kerrebejin et al. GWAS of Fibromyalgia[@kerrebijn2026genetic] colocalized with proteomics GWAS from the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md)[@sun2023plasma].


Locus/protein pairs with greater than 80 percent probability of colocalization are shown below[^filtering]


{{
coloc_results_table("assets/base_asset_store/gwas/fibromyalgia/kerrebijn_et_al/analysis/kerrebijn_fibro_ppp_coloc_full_ppp_coloc_display.parquet","fibro_ppp_coloc")
}}


















[^filtering]: I used a filtering algorithm to select which protein/locus pairs to test.  Consequently, some pairs that looked a priori unlikely to produce useful signal were not tested.