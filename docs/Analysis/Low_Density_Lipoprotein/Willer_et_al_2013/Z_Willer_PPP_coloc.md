---
tags:
- ABF colocalization
---
# PPP Coloc

For further insights into the mechanisms determining LDL levels, I ran [ABF Colocalization](../../../Bioinformatics_Concepts/Colocalization.md)[@giambartolomei2014bayesian] to search for loci at which Willer et al.'s GWAS of LDL levels[@global2013discovery] colocalized with proteomics GWAS from the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md)[@sun2023plasma]. Because PPP GWAS files are large and slow to download, I limited the scan to 30 protein/locus pairs[^scan_note].


The results are below.


{{
coloc_results_table("docs/_figs/willer_ldl_ppp_coloc_display.parquet","willer_ppp_coloc",)
}}


- PCSK9 is an expected and highly scientifically-supported result[^podcast_note].






[^scan_note]: I selected these pairs by prioritizing proteins/loci for which both LDL and the protein had a low p value at the locus.


[^podcast_note]: For example, see [this podcast episode](https://podcasts.apple.com/us/podcast/should-everyone-be-taking-statins/id1819489230?i=1000751897974) which includes discussions of PCSK9 as breakthrough LDL medications.
                           


