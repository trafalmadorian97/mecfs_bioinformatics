# S-LDSC

I applied stratified linkage disequilibrium score regression (S-LDSC)[@finucane2018heritability] to summary statistics from DECODE's meta-GWAS of seropositive rheumatoid arthritis (RA)[@saevarsdottir2022multiomics].


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

The plot and expandable table below show the results of S-LDSC applied to the seropositive RA GWAS using the GTEx and Franke lab gene expression reference datasets. In the plot, the x-axis corresponds to cell type, while the y-axis corresponds to $-\log_{10}(p)$.  Points are colored according to tissue category.  Large points indicate cell/tissue types deemed significant by the Benjamini-Hochberg procedure at an FDR of 0.01[@benjamini1995controlling]. 

{{ plotly_embed("docs/_figs/decode_ra_seropositive_multi_tissue_gene_expression_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="seropositive-ra-sldsc-gene-expression") }}


{{ markdown_table("docs/_figs/decode_ra_seropositive_multi_tissue_gene_expression_s_ldsc_cell_analysis_md_table.mdx", title="GTEx and Franke lab tissue expression — full results", collapse_threshold=0) }}

As would be expected for an immune trait, the main significant cells/tissue are immune-related. There is one significant gut tissue type. It is unclear how to interpret this, but it may reflect the presence of immune cells in the gut.


### Roadmap/ ENCODE Chromatin data

I next applied S-LDSC using the Roadmap/ ENCODE epigenetic reference dataset.  The results are in the plot and expandable table below:

{{ plotly_embed("docs/_figs/decode_ra_seropositive_multi_tissue_chromatin_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="decode-ra-sldsc-chromatin") }}

{{ markdown_table("docs/_figs/decode_ra_seropositive_multi_tissue_chromatin_s_ldsc_cell_analysis_md_table.mdx", title="Roadmap/ ENCPDE chromatin — full results", collapse_threshold=0) }}


Epigenic reference data reinforces the point made by the GTEx and Franke lab reference data: the key cell types pertain to the immune system.


### ImmGen data

Next, I applied S-LDSC using the ImmGen project reference dataset.

{{ plotly_embed("docs/_figs/decode_ra_seropositive_immgen_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="decode-ra-sldsc-immgen") }}

{{ markdown_table("docs/_figs/decode_ra_seropositive_immgen_s_ldsc_cell_analysis_md_table.mdx", title="Immgen — full results", collapse_threshold=0) }}

Consistent with rheumatoid arthritis's known status as a disease of the adaptive immune system, we see that T and B cells are significant.  The significance of natural killer cells is interesting, but may be explained by their transcriptional similarity to T cells.

### Corces et al. ATAC-seq data

The results of applying S-LDSC using the epigenetic reference data from the Corces et al. ATAC-seq analysis of blood cells are shown in the expandable table below.

{{ markdown_table("docs/_figs/decode_ra_seropositive_corces_atac_s_ldsc_cell_analysis_md_table.mdx", title="Corces et al. — full results", collapse_threshold=0) }}

Again, T cells, B cells, and natural killer cells are significant, consistent with an adaptive-immune pathobiology.


### Cahoy and GTEx-Brain data

The next two reference datasets pertain to the nervous system.



{{ markdown_table("docs/_figs/decode_ra_seropositive_cahoy_cns_s_ldsc_cell_analysis_md_table.mdx", title="Cahoy — full results", collapse_threshold=0) }}



{{ markdown_table("docs/_figs/decode_ra_seropositive_gtex_brain_s_ldsc_cell_analysis_md_table.mdx", title="GTEx brain — full results", collapse_threshold=0) }}

Note of the cells or tissue types in these datasets are significant.  This is consistent with RA being a non-neural disease.


