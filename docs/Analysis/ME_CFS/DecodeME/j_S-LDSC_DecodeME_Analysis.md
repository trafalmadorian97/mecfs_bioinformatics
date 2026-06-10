---
hide:
- toc
---
# S-LDSC

[Stratified Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) (S-LDSC)[@finucane2018heritability] was applied to summary statistics from GWAS-1 of [DecodeME](../../../Data_Sources/DecodeME.md)[@genetics2025initial].

## Reference Data Sources

I used the standard reference datasets prepared by the authors of the S-LDSC method.

- [The GTEx Project](../../../Data_Sources/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../../Data_Sources/Roadmap.md)
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../../Data_Sources/Corces_et_al.md).
- The [ImmGen](../../../Data_Sources/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset


## Results

### GTEx and Franke lab tissue expression data
The plot and table below show the results of the application of S-LDSC to DecodeME using the GTEx and Franke lab gene expression datasets. In the plot, the x-axis corresponds to cell type, while the y-axis corresponds to $-\log_{10}(p)$.  Points are colored according to broad tissue category.  Large points correspond to cell/tissue types deemed significant by an application of the Benjamini-Hochberg procedure at an FDR of 0.01[@benjamini1995controlling].  The table lists every cell/tissue type in the reference panel, sorted by p value.


{{ plotly_embed("docs/_figs/decode_me_gwas_1_multi_tissue_gene_expression_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="decode-me-sldsc-gene-expression") }}

{{ markdown_table("docs/_figs/decode_me_gwas_1_multi_tissue_gene_expression_s_ldsc_cell_analysis_md_table.mdx", title="GTEx and Franke lab tissue expression — full results") }}


As we saw in our earlier [MAGMA analysis](g_MAGMA_DecodeME_Analysis.md) using the GTEx dataset, the significant tissues are all CNS-related.


### Roadmap Chromatin data

I next applied S-LDSC using the reference dataset derived from the Roadmap epigenetic project.  The results are in the plot and table below:

{{ plotly_embed("docs/_figs/decode_me_gwas_1_multi_tissue_chromatin_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="decode-me-sldsc-chromatin") }}

{{ markdown_table("docs/_figs/decode_me_gwas_1_multi_tissue_chromatin_s_ldsc_cell_analysis_md_table.mdx", title="Roadmap chromatin — full results") }}

Again, the strongest and most significant associations are all with CNS cell-types.

### ImmGen data

Next, I applied S-LDSC using reference data from the ImmGen project.

There were no significant cell types.

The full set of cell types, sorted by p value, is shown in the table below:

{{ markdown_table("docs/_figs/decode_me_gwas_1_immgen_s_ldsc_cell_analysis_md_table.mdx", title="ImmGen — full results") }}


### Corces et al. ATAC-seq data

The results of applying S-LDSC using the epigenetic reference data from Corces et al. ATAC-seq analysis of blood cells are shown below.  There are no significant cell types:

{{ markdown_table("docs/_figs/decode_me_gwas_1_corces_atac_s_ldsc_cell_analysis_md_table.mdx", title="Corces et al. ATAC-seq — full results") }}


### Cahoy and GTEx-Brain data

The next two reference datasets pertain to the nervous system.

Surprisingly, when we analyze the DecodeME results using the Cahoy dataset, the neuron cell type is not even close to being significant.  This is discordant with some of the results above, in which many CNS-related cell and tissue types were marked as significant.  Moreover, the oligodendrocyte cell type is closer to being significant than the neuron cell type.

{{ markdown_table("docs/_figs/decode_me_gwas_1_cahoy_cns_s_ldsc_cell_analysis_md_table.mdx", title="Cahoy — full results") }}


When we apply the S-LDSC using the GTEx brain dataset, we find the the cortex tissue type is significant:

{{ markdown_table("docs/_figs/decode_me_gwas_1_gtex_brain_s_ldsc_cell_analysis_md_table.mdx", title="GTEx brain — full results") }}


## How to Reproduce This

To reproduce, run the {{ api_link("DecodeME Analysis Script", "mecfs_bio.analysis.decode_me_initial_analysis") }}.
