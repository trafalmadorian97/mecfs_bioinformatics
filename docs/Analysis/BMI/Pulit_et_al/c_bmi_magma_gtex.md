---
tags:
- MAGMA
---
# MAGMA GTEx


I applied [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to the Pulit et al. meta-GWAS of BMI[@pulit2019meta]. 



## MAGMA Gene Analysis

I used MAGMA's SNP-wise mean model to run gene-level analysis.


In this step:

- Data from the 1000 genomes projects was downloaded from the [MAGMA website](https://cncr.nl/research/magma/) and used as a linkage disequilibrium reference.
- Magma's default proximity-based rules were used to assign SNPs to genes.


The plot below shows the results of MAGMA gene-level analysis

{{ plotly_embed("docs/_figs/pulit_et_al_2018_bmi_magma_gene_manhattan_plot.html", id="magma-bmi-gene-manhattan-gene-manhattan-window") }}



## MAGMA Gene Property Analysis


I next applied [MAGMA's](../../../Bioinformatics_Concepts/MAGMA_Overview.md) gene property analysis module to the BMI GWAS.  This step combines the gene analysis results above with tissue-specific RNA expression values from the [GTEx project](../../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex]. The resulting  per-tissue p values are plotted below:


{{ plotly_embed("docs/_figs/pulit_et_al_2018_bmi_magma_bar_plot/magma_gene_set_plot.html", id="bmi-magma-gtex-tissue-bar-plot") }}


The strong enrichment of brain tissue is consistent the known importance of neural pathways in obesity[^gene_eating].

[^gene_eating]: See Yeo 2018[@yeo2018gene] for a readable popular science account of the genetics of eating.