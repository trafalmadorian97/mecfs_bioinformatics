---
hide:
- toc
---
# S-LDSC

I applied [Stratified Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) (S-LDSC)[@finucane2018heritability] to summary statistics from [DecodeME](../../../Data_Sources/DecodeME.md)[@genetics2025initial] GWAS-1.

## Reference Data Sources

I used the standard reference datasets prepared by the authors of the S-LDSC method.

- [The GTEx Project](../../../Data_Sources/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../../Data_Sources/Roadmap.md) and the ENCODE project.
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../../Data_Sources/Corces_et_al.md).
- The [ImmGen](../../../Data_Sources/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset


## Results

### GTEx and Franke lab tissue expression data
The plot and expandable table below show the results of S-LDSC applied to DecodeME using the GTEx and Franke lab gene expression reference datasets. In the plot, the x-axis corresponds to cell type, while the y-axis corresponds to $-\log_{10}(p)$.  Points are colored according to tissue category.  Large points indicate cell/tissue types deemed significant by the Benjamini-Hochberg procedure at an FDR of 0.01[@benjamini1995controlling].  


{{ plotly_embed("docs/_figs/decode_me_gwas_1_multi_tissue_gene_expression_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="decode-me-sldsc-gene-expression") }}

{{ markdown_table("docs/_figs/decode_me_gwas_1_multi_tissue_gene_expression_s_ldsc_cell_analysis_md_table.mdx", title="GTEx and Franke lab tissue expression — full results", collapse_threshold=0) }}


As we saw in [GTEx MAGMA analysis](g_MAGMA_DecodeME_Analysis.md), the significant tissues are all CNS-related.


### Roadmap/ ENCODE Chromatin data

I next applied S-LDSC using the Roadmap epigenetic reference dataset.  The results are in the plot and expandable table below:

{{ plotly_embed("docs/_figs/decode_me_gwas_1_multi_tissue_chromatin_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="decode-me-sldsc-chromatin") }}

{{ markdown_table("docs/_figs/decode_me_gwas_1_multi_tissue_chromatin_s_ldsc_cell_analysis_md_table.mdx", title="Roadmap chromatin — full results", collapse_threshold=0) }}

Again, the strongest associations are with CNS cell-types.

### ImmGen data

Next, I applied S-LDSC using the ImmGen project reference dataset.

There were no significant cell types.

The full results are in the expandable table below.

{{ markdown_table("docs/_figs/decode_me_gwas_1_immgen_s_ldsc_cell_analysis_md_table.mdx", title="ImmGen — full results", collapse_threshold=0) }}


### Corces et al. ATAC-seq data

The results of applying S-LDSC using the epigenetic reference data from Corces et al. ATAC-seq analysis of blood cells are shown in the expandable table below.  There are no significant cell types.

{{ markdown_table("docs/_figs/decode_me_gwas_1_corces_atac_s_ldsc_cell_analysis_md_table.mdx", title="Corces et al. ATAC-seq — full results", collapse_threshold=0) }}


### Cahoy and GTEx-Brain data

The next two reference datasets pertain to the nervous system.

Surprisingly, when we analyze the DecodeME results using the Cahoy dataset, the neuron cell type is not even close to being significant.  This is discordant with some of the results above, in which many CNS-related cell and tissue types were marked as significant.  Moreover, the oligodendrocyte cell type is closer to being significant than the neuron cell type.

{{ markdown_table("docs/_figs/decode_me_gwas_1_cahoy_cns_s_ldsc_cell_analysis_md_table.mdx", title="Cahoy — full results", collapse_threshold=0) }}


When we use the GTEx brain reference dataset, the cortex tissue type is significant. See the expandable table.

{{ markdown_table("docs/_figs/decode_me_gwas_1_gtex_brain_s_ldsc_cell_analysis_md_table.mdx", title="GTEx brain — full results", collapse_threshold=0) }}


## How to Reproduce This

To reproduce, run the {{ api_link("DecodeME Analysis Script", "mecfs_bio.analysis.decode_me_initial_analysis") }}.
