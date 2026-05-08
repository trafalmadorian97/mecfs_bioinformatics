---
hide:
- toc
---
# MAGMA



I applied the [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] pipeline to summary statistics from a GWAS of C-Reactive Protein (CRP) from Said at el.[@said2022genetic].


In the gene set analysis step, I incorporated [tissue-specific RNAseq data from GTEx](../../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex] to try to link genes associated with LDL to specific tissues.

The results are shown below:

{{ plotly_embed("docs/_figs/said_et_al_2022_crp_eur_magma_bar_plot/magma_gene_set_plot.html", id="said-crp-gtex-magma", caption="Result of MAGMA GTEx applied to Said et al.'s CRP GWAS. The x axis corresponds to tissue type, while the y axis measures MAGMA significance level.") }}


All the significant tissue are liver-related.  This is consistent with the liver's role as a primary site of acute-phase-reactant synthesis.
