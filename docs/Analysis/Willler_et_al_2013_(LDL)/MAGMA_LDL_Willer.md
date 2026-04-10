---
hide:
- toc
---
# MAGMA 

I applied the [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] pipeline to summary statistics from a GWAS of LDL from WIller at el.[@global2013discovery].

 In the gene set analysis step, I incorporated [tissue-specific RNAseq data from GTEx](../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex] to try to link genes associated with LDL to specific tissues.

The results are shown below:



{{ plotly_embed("docs/_figs/willer_et_al_2023_ldl_eur_magma_bar_plot/magma_gene_set_plot.html", id="willer-ldl-gtex-magma", caption="Result of MAGMA GTEx applied to Willer et al.'s LDL GWAS. The x axis corresponds to tissue type, while the y axis measures MAGMA significance level.") }}


The 4 most-significant tissues found by MAGMA are all in liver.  This is very consistent with the known biology of LDL[^biology_ldl_note] which is both synthesized and disposed of by liver cells.



[^biology_ldl_note]: See for instance Chapters 4 and 5 of Steinberg's book[@steinberg2011cholesterol].