---
tags:
- ABF colocalization
---
# PPP Coloc
For further insights into the mechanisms of seropositive rheumatoid arthritis, I ran [ABF Colocalization](../../../Bioinformatics_Concepts/Colocalization.md)[@giambartolomei2014bayesian] to search for loci at which the DECODE meta-GWAS of seropositive RA[@saevarsdottir2022multiomics] colocalized with proteomics GWAS from the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md)[@sun2023plasma]. Because PPP GWAS files are large and slow to download, I limited the scan to 30 protein/locus pairs


The results are below.


{{
coloc_results_table("docs/_figs/decode_ra_seropositive_ppp_coloc_display.parquet","ra_ppp_coloc")
}}


There is a single colocalized protein: CCL19.  Thus we find that in seropositive rheumatoid arthritis, both the global polygenic signal and the locus-specific signal point to shared genetic architecture with blood levels of CCL19.