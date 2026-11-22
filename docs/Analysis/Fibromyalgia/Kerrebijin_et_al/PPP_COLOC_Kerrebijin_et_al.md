# PPP Coloc
For further insights into the molecular mechanisms of Fibromyalgia, I ran [ABF Colocalization](../../../Bioinformatics_Concepts/Colocalization.md)[@giambartolomei2014bayesian] to search for loci at which the Kerrebejin et al. GWAS of Fibromyalgia[@kerrebijn2026genetic] colocalized with proteomics GWAS from the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md)[@sun2023plasma].


Locus/protein pairs with greater than 80 percent probability of colocalization are shown below


{{
coloc_results_table("docs/_figs/kerrebijn_fibro_ppp_coloc_full_ppp_coloc_display.parquet","fibro_ppp_coloc")
}}



- The only cis colocalization is with [SUOX](https://www.genecards.org/card/SUOX) at the locus **chr12:55491020-56491020**.  At this locus, fibromyalgia also colocalized with numerous other proteins in trans. Notably, this included [IL12B](https://www.genecards.org/card/IL12B), a subunit of the cytokine interleukin 12.













