## MAGMA DecodeME Analysis
As an initial step, I applied [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md) to [DecodeME](../../Data_Sources/MECFS/DecodeME.md), partially reproducing analysis from the DecodeME preprint.  


## MAGMA Gene Analysis

I applied MAGMA's SNP-wise mean model to the summary statistics from DecodeME's GWAS 1. 


In this step:

- Data from the 1000 genomes projects was downloaded from the [MAGMA website](https://cncr.nl/research/magma/) and used as a linkage disequilibrium reference.
- Data from the [SNP151 database](https://hgw2.soe.ucsc.edu/cgi-bin/hgTables?hgsid=2912494930_cRufLdpdc1ynRc2sCM3g1WGAWAgH&hgta_doSchemaDb=hg19&hgta_doSchemaTable=snp151Flagged
  ) was used to assign RSIDs to SNPs.
- Magma's default proximity-based rules were used to assign SNPs to fenes.

MAGMA produces a table of genes, effect sizes, and p values.  Filtering these genes via the Benjamini-Hochberg procedure[@benjamini1995controlling] at a false discovery rate of 0.01 produces:

 | GENE            |   CHR |     START |      STOP |   NSNPS |   NPARAM |      N |   ZSTAT |          P |   _Corrected P Value_ |
|:----------------|------:|----------:|----------:|--------:|---------:|-------:|--------:|-----------:|----------------------:|
| ENSG00000033122 |     1 |  70034081 |  70617628 |    1609 |       68 | 275488 |  5.9087 | 1.724e-09  |           3.22647e-05 |
| ENSG00000124214 |    20 |  47729878 |  47804904 |     152 |       16 | 275488 |  5.6715 | 7.0767e-09 |           6.62202e-05 |
| ENSG00000124207 |    20 |  47662849 |  47713489 |     116 |       11 | 275488 |  5.4047 | 3.2463e-08 |           0.000202515 |
| ENSG00000135090 |    12 | 118587606 | 118810750 |     445 |       24 | 275488 |  5.0416 | 2.3082e-07 |           0.00104115  |
| ENSG00000124198 |    20 |  47538427 |  47653230 |     262 |       17 | 275488 |  4.9822 | 3.1435e-07 |           0.00104115  |
| ENSG00000117593 |     1 | 173793641 | 173827684 |      37 |        7 | 275488 |  4.9706 | 3.3379e-07 |           0.00104115  |
| ENSG00000185278 |     1 | 173837220 | 173872687 |     105 |       13 | 275488 |  4.9122 | 4.503e-07  |           0.00120391  |
| ENSG00000250091 |    12 | 124410971 | 124419531 |      14 |        3 | 275488 |  4.8395 | 6.5098e-07 |           0.00152289  |
| ENSG00000179195 |    12 | 124456392 | 124499986 |      77 |        9 | 275488 |  4.7302 | 1.1217e-06 |           0.00219969  |
| ENSG00000119242 |    12 | 124403207 | 124457378 |     116 |        9 | 275488 |  4.7152 | 1.2073e-06 |           0.00219969  |
| ENSG00000197935 |     6 |  28962562 |  28973093 |      30 |        8 | 275488 |  4.7013 | 1.2929e-06 |           0.00219969  |
| ENSG00000158406 |     6 |  26281283 |  26285762 |      13 |        4 | 275488 |  4.6324 | 1.8071e-06 |           0.00281832  |
| ENSG00000198216 |     1 | 181382238 | 181777219 |    1030 |       79 | 275488 |  4.591  | 2.2062e-06 |           0.00299373  |
| ENSG00000187323 |    18 |  49866542 |  51057784 |    4538 |      112 | 275488 |  4.5878 | 2.2395e-06 |           0.00299373  |
| ENSG00000111707 |    12 | 118814185 | 118855840 |     114 |       11 | 275488 |  4.5649 | 2.4987e-06 |           0.00311754  |
| ENSG00000197653 |    12 | 124247042 | 124420753 |     458 |       42 | 275488 |  4.4981 | 3.4281e-06 |           0.00400981  |
| ENSG00000089220 |    12 | 118573663 | 118583389 |      29 |        8 | 275488 |  4.4695 | 3.9196e-06 |           0.00431502  |
| ENSG00000188730 |     7 |  49813257 |  49961546 |     320 |       34 | 275488 |  4.3858 | 5.7778e-06 |           0.00600731  |
| ENSG00000028116 |     2 |  58134786 |  58387055 |     503 |       26 | 275488 |  4.33   | 7.4569e-06 |           0.00734505  |
| ENSG00000117601 |     1 | 173872947 | 173886516 |      26 |        5 | 275488 |  4.2942 | 8.7656e-06 |           0.00820241  |

Note that the genomic coordinates in the table above refer to genome build 37.

[//]: # (To reproduce this analysis use the build system to materialize the asset declared in [this file]&#40;https://github.com/trafalmadorian97/mecfs_bioinformatics/blob/ba3ee982c59977358ee0d40708347d43a009f190/src_new/assets/gwas/me_cfs/decode_me/processed_gwas_data/magma/decode_me_gwas_1_build_37_magma_ensembl_gene_analysis.py&#41;.)

## MAGMA Gene Property Analysis

I next applied [MAGMA's](../../Bioinformatics_Concepts/MAGMA_Overview.md) gene property analysis module to [DecodeME](../../Data_Sources/MECFS/DecodeME.md).  This step combined the gene analysis results above with tissue-specific RNA expression values from the [GTEx project](../../Data_Sources/GTEx_project/GTEx_RNAseq_Data.md).  The aim was to identify tissues enriched for genes associated with ME/CFS.  The results are shown in the bar plot below:

![bar_plot_decodeme_tissues](https://github.com/user-attachments/assets/519629d4-c60a-434b-9ceb-809c2878cbe5)
In this plot, the y axis corresponds to negative log p values, the x axis corresponds to tissue type (only tissues with the smallest p values are shown), and bars are colored according to whether their p-value meets Bonferroni-corrected significance threshold,



These results unambiguously point to the nervous system as a major site of ME/CFS gene activity.

[//]: # (To reproduce this analysis, use the build system to materialize the asset in [this file]&#40;https://github.com/trafalmadorian97/mecfs_bioinformatics/blob/ba3ee982c59977358ee0d40708347d43a009f190/src_new/assets/gwas/me_cfs/decode_me/analysis_results/magma/magma_specific_tissue_bar_plot.py&#41;.)


## Reproducing
To reproduce this analysis, run the [DecodeME Initial Analysis Script][mecfs_bio.analysis.decode_me_initial_analysis].

## Follow-Up Questions
1. Do other approaches to identify significant tissues from GWAS-summary statistics produce concordant or discordant results?
2. How reliable is the MAGMA gene-set-analysis approach to identifying significant tissues?  In other words: for diseases with well-understood pathological processes, does it produce results consistent with these processes?