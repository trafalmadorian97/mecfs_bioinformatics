## MAGMA DecodeME Analysis
As an initial step, I applied [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md) to [DecodeME](../../Data_Sources/MECFS/DecodeME.md), partially reproducing analysis from the DecodeME preprint.  


## MAGMA Gene Analysis

I applied MAGMA's SNP-wise mean model to the summary statistics from DecodeME's GWAS 1. 


In this step:

- Data from the 1000 genomes projects was downloaded from the [MAGMA website](https://cncr.nl/research/magma/) and used as a linkage disequilibrium reference.
- Data from the [SNP151 database](https://hgw2.soe.ucsc.edu/cgi-bin/hgTables?hgsid=2912494930_cRufLdpdc1ynRc2sCM3g1WGAWAgH&hgta_doSchemaDb=hg19&hgta_doSchemaTable=snp151Flagged
  ) was used to assign RSIDs to SNPs.
- Magma's default proximity-based rules were used to assign SNPs to fenes.

MAGMA produces a table of genes, effect sizes, and p values.  Filtering these genes via the Benjamini-Hochberg procedure[@benjamini1995controlling] at a false discovery rate of 0.01, and joining with a database of gene descriptions from [Ensembl Biomart](https://useast.ensembl.org/info/data/biomart/index.html) produces the following table:

| GENE            | Gene name   |   CHR |          P | Gene description                                                                      |
|:----------------|:------------|------:|-----------:|:--------------------------------------------------------------------------------------|
| ENSG00000033122 | LRRC7       |     1 | 1.724e-09  | leucine rich repeat containing 7 [Source:HGNC Symbol;Acc:HGNC:18531]                  |
| ENSG00000124214 | STAU1       |    20 | 7.0767e-09 | staufen double-stranded RNA binding protein 1 [Source:HGNC Symbol;Acc:HGNC:11370]     |
| ENSG00000124207 | CSE1L       |    20 | 3.2463e-08 | chromosome segregation 1 like [Source:HGNC Symbol;Acc:HGNC:2431]                      |
| ENSG00000135090 | TAOK3       |    12 | 2.3082e-07 | TAO kinase 3 [Source:HGNC Symbol;Acc:HGNC:18133]                                      |
| ENSG00000124198 | ARFGEF2     |    20 | 3.1435e-07 | ARF guanine nucleotide exchange factor 2 [Source:HGNC Symbol;Acc:HGNC:15853]          |
| ENSG00000117593 | DARS2       |     1 | 3.3379e-07 | aspartyl-tRNA synthetase 2, mitochondrial [Source:HGNC Symbol;Acc:HGNC:25538]         |
| ENSG00000185278 | ZBTB37      |     1 | 4.503e-07  | zinc finger and BTB domain containing 37 [Source:HGNC Symbol;Acc:HGNC:28365]          |
| ENSG00000250091 | DNAH10OS    |    12 | 6.5098e-07 | dynein axonemal heavy chain 10 opposite strand [Source:HGNC Symbol;Acc:HGNC:37121]    |
| ENSG00000179195 | ZNF664      |    12 | 1.1217e-06 | zinc finger protein 664 [Source:HGNC Symbol;Acc:HGNC:25406]                           |
| ENSG00000119242 | CCDC92      |    12 | 1.2073e-06 | coiled-coil domain containing 92 [Source:HGNC Symbol;Acc:HGNC:29563]                  |
| ENSG00000197935 | ZNF311      |     6 | 1.2929e-06 | zinc finger protein 311 [Source:HGNC Symbol;Acc:HGNC:13847]                           |
| ENSG00000158406 | H4C8        |     6 | 1.8071e-06 | H4 clustered histone 8 [Source:HGNC Symbol;Acc:HGNC:4788]                             |
| ENSG00000198216 | CACNA1E     |     1 | 2.2062e-06 | calcium voltage-gated channel subunit alpha1 E [Source:HGNC Symbol;Acc:HGNC:1392]     |
| ENSG00000187323 | DCC         |    18 | 2.2395e-06 | DCC netrin 1 receptor [Source:HGNC Symbol;Acc:HGNC:2701]                              |
| ENSG00000111707 | SUDS3       |    12 | 2.4987e-06 | SDS3 homolog, SIN3A corepressor complex component [Source:HGNC Symbol;Acc:HGNC:29545] |
| ENSG00000197653 | DNAH10      |    12 | 3.4281e-06 | dynein axonemal heavy chain 10 [Source:HGNC Symbol;Acc:HGNC:2941]                     |
| ENSG00000089220 | PEBP1       |    12 | 3.9196e-06 | phosphatidylethanolamine binding protein 1 [Source:HGNC Symbol;Acc:HGNC:8630]         |
| ENSG00000188730 | VWC2        |     7 | 5.7778e-06 | von Willebrand factor C domain containing 2 [Source:HGNC Symbol;Acc:HGNC:30200]       |
| ENSG00000028116 | VRK2        |     2 | 7.4569e-06 | VRK serine/threonine kinase 2 [Source:HGNC Symbol;Acc:HGNC:12719]                     |
| ENSG00000117601 | SERPINC1    |     1 | 8.7656e-06 | serpin family C member 1 [Source:HGNC Symbol;Acc:HGNC:775]                            |


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