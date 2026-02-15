---
hide:
  - navigation
  - toc
---
# Chr20 Locus

In an attempt to narrow down the DecodeME[@genetics2025initial] GWAS-1 signal, I [fine-mapped](../../../Bioinformatics_Concepts/Fine_Mapping.md) the GWAS-1 hit on chromosome 20 using SUSIE[@wang2020simple].


As a linkage disequilibrium reference, I used a UK Biobank LD matrix hosted on [AWS Open Data](https://registry.opendata.aws/ukbb-ld/).  Because this LD reference uses GRCh37 coordinates, I used liftover to translate the DecodeME GWAS-1 summary statistics to GRCh37.

## Results with $L=10$

The results of running SUSIE on the chromosome 20 DecodeME GWAS-1 locus are plotted below:


![decode_me_susie_chr20](../../../_figs/decode_me_chr20_46000001_49000001_susie_stackplot.png)


SUSIE reports that it has found two credible sets, one of which is diffuse, and one of which is highly concentrated.