---
hide:
  - toc
---
# MAGMA

I applied the [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] pipeline to summary statistics from a GWAS of triglycerides from Willer at el.[@global2013discovery].


In the gene set analysis step, I incorporated [tissue-specific RNAseq data from GTEx](../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex] to try to link genes associated with LDL to specific tissues.

The results are shown below:

{{ plotly_embed("../../../_figs/willer_et_al_2023_tg_eur_magma_bar_plot/magma_gene_set_plot.html", id="willer-tg-gtex-magma", caption="Result of MAGMA GTEx applied to Willer et al.'s triglyceride GWAS. The x axis corresponds to tissue type, while the y axis measures MAGMA significance level.") }}


All the significant tissue are liver-related.  This is consistent with the liver's role as a central site of lipid metabolism.