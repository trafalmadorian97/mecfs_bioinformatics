---
hide:
- navigation
- toc
---
# Chr1 Locus

In an attempt to narrow down the DecodeME[@genetics2025initial] GWAS-1 signal, I [fine-mapped](../../../../Bioinformatics_Concepts/Fine_Mapping.md) the GWAS-1 hit on chromosome 1 using SUSIE[@wang2020simple].

As a linkage disequilibrium reference, I used a UK Biobank LD matrix hosted on [AWS Open Data](https://registry.opendata.aws/ukbb-ld/).  Because this LD reference uses GRCh37 coordinates, I used liftover to translate the DecodeME GWAS-1 summary statistics to GRCh37.

The results of running SUSIE on the chromosome 1 DecodeME GWAS-1 locus are plotted below:

[//]: # (![decodeme_ch1_stackplot]&#40;../../../_figs/stackplot_decode_me_gwas_1_37_susie_finemap_chr1_173000001_locus.png&#41;)

- The top panel is a heatmap in pixel (i,j) is colored according to the  correlation of between variants i and j.  The heatmap thus reveals the local linkage disequilibrium (LD) structure in the vicinity of the GWAS hit.

- The second panel shows a local Manhattan plot.

- The third panel shows SUSIE posterior inclusion probability (PIP).

- The bottom panel shows genes in the region of the GWAS hit. 

In this instance, SUSIE has found a single credible set containing many variants, each with small PIP.  This non-specific finding may be due to the GWAS hit lying within a large LD block, making it difficult for SUSIE to narrow down a specific variant.  It may be worth trying other fine-mapping algorithms in this locus.