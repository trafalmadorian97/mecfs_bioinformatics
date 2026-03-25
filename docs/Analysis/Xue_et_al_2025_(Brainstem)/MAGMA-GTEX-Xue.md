---
hide:
- toc
---
# MAGMA GTEx

I applied [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to the whole brainstem volume GWAS of Xue et al.[@xue2025genetic] using bulk RNAseq data from [GTEx](../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex] as a reference. The gene-property analysis results are plotted below.

{{ plotly_embed("../../../_figs/xue_whole_brainstem_magma_bar_plot/magma_gene_set_plot.html", id="xue-whole-brainstem-gtex-magma", caption="Result of MAGMA GTEx applied to Xue et al.'s whole brainstem GWAS. The x axis corresponds to tissue type, while the y axis measures MAGMA significance level.") }}

As would be expected from a GWAS of the volume of a region of the brain, all 4 significant tissue types relate to the nervous system.  These consist of three central nervous system cell types: Brain_Cerebellum, Brain_Cerebellar_Hemisphere, and Brain_Spinal_cord_cervical_c-1; as well as one peripheral nervous system cell type: Nerve_Tibial.  It is unclear why these specific nervous system cell types are selected by MAGMA over other nerovus-system related cell types.