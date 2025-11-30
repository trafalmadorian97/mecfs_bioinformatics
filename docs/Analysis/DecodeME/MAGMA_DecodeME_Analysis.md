## MAGMA DecodeME Analysis
As an initial step, I applied [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md) to [DecodeME](../../Data_Sources/MECFS/DecodeME.md), partially reproducing analysis from the DecodeME preprint.  


## MAGMA Gene Analysis

I applied MAGMA's SNP-wise mean model to the summary statistics from DecodeME's GWAS 1. 


In this step:

- Data from the 1000 genomes projects was downloaded from the [MAGMA website](https://cncr.nl/research/magma/) and used as a linkage disequilibrium reference.
- Data from the [SNP151 database](https://hgw2.soe.ucsc.edu/cgi-bin/hgTables?hgsid=2912494930_cRufLdpdc1ynRc2sCM3g1WGAWAgH&hgta_doSchemaDb=hg19&hgta_doSchemaTable=snp151Flagged
  ) was used to assign RSIDs to SNPs.
- Magma's default proximity-based rules were used to assign SNPs to fenes.

MAGMA produces a table of genes, effect sizes, and p values of the form:
```
GENE             CHR      START       STOP  NSNPS  NPARAM       N        ZSTAT            P
ENSG00000269831    1     738532     739137      1       1  275488      0.53146      0.29755
ENSG00000187634    1     860260     879955     63      15  275488      0.93486      0.17493
ENSG00000268179    1     861264     866445     24       7  275488       1.5318     0.062784
ENSG00000188976    1     879584     894689     44       7  275488      0.47458      0.31754
ENSG00000187961    1     895967     901095     16       4  275488      0.13407      0.44667
ENSG00000187583    1     901877     911245     33      14  275488       1.0717      0.14193
ENSG00000187642    1     910579     917497     22       5  275488       1.1909      0.11685
ENSG00000188290    1     934342     935552      5       2  275488        1.243      0.10693
ENSG00000187608    1     948803     949920      5       1  275488    -0.057976      0.52312
ENSG00000188157    1     955503     991496     94      12  275488      -2.1128      0.98269
ENSG00000237330    1    1006346    1009687      6       2  275488     -0.77358      0.78041
...
```




To reproduce this analysis use the build system to materialize the asset declared in [this file](https://github.com/trafalmadorian97/mecfs_bioinformatics/blob/ba3ee982c59977358ee0d40708347d43a009f190/src_new/assets/gwas/me_cfs/decode_me/processed_gwas_data/magma/decode_me_gwas_1_build_37_magma_ensembl_gene_analysis.py).

## MAGMA Gene Property Analysis

I next applied [MAGMA's](../../Bioinformatics_Concepts/MAGMA_Overview.md) gene property analysis module to [DecodeME](../../Data_Sources/MECFS/DecodeME.md).  This step combined the gene analysis results above with tissue-specific RNA expression values from the [GTEx project](../../Data_Sources/GTEx_project/GTEx_RNAseq_Data.md).  The aim was to identify tissues enriched for genes associated with ME/CFS.  The results are shown in the bar plot below:

![bar_plot_decodeme_tissues](https://github.com/user-attachments/assets/519629d4-c60a-434b-9ceb-809c2878cbe5)
In this plot, the y axis corresponds to negative log p values, the x axis corresponds to tissue type (only tissues with the smallest p values are shown), and bars are colored according to whether their p-value meets Bonferroni-corrected significance threshold,



These results unambiguously point to the nervous system as a major site of ME/CFS gene activity.

To reproduce this analysis, use the build system to materialize the asset in [this file](https://github.com/trafalmadorian97/mecfs_bioinformatics/blob/ba3ee982c59977358ee0d40708347d43a009f190/src_new/assets/gwas/me_cfs/decode_me/analysis_results/magma/magma_specific_tissue_bar_plot.py).


## Follow-Up Questions
1. Do other approaches to identify significant tissues from GWAS-summary statistics produce concordant or discordant results?
2. How reliable is the MAGMA gene-set-analysis approach to identifying significant tissues?  In other words: for diseases with well-understood pathological processes, does it produce results consistent with these processes?