---
tags:
- S-LDSC
---
# S-LDSC

I applied stratified linkage disequilibrium score regression (S-LDSC)[@finucane2018heritability] to summary statistics from Pulet et al. meta-GWAS of body mass index (BMI)[@pulit2019meta].



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



The plot and table below show the results of S-LDSC applied to the seropositive RA GWAS using the GTEx and Franke lab gene expression reference datasets. In the plot, the x-axis corresponds to cell type, while the y-axis corresponds to $-\log_{10}(p)$.  Points are colored according to tissue category.  Large points indicate cell/tissue types deemed significant by the Benjamini-Hochberg procedure at an FDR of 0.01[@benjamini1995controlling].

{{ plotly_embed("docs/_figs/pulit_et_al_2018_bmi_multi_tissue_gene_expression_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="bmi-sldsc-gene-expression") }}


{{
data_table("docs/_figs/pulit_et_al_2018_bmi_multi_tissue_gene_expression_s_ldsc_cell_analysis_parquet_table.parquet", id="bmi-gene-expression-table")
}}

Consistent with the known biology of obesity and our previous [MAGMA results](c_bmi_magma_gtex.md), all significant associations are in the central nervous system category.


### Roadmap/ ENCODE Chromatin data

I next applied S-LDSC using the Roadmap/ ENCODE epigenetic reference dataset.  The results are in the plot and expandable table below:

{{ plotly_embed("docs/_figs/pulit_et_al_2018_bmi_multi_tissue_chromatin_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="decode-bmi-sldsc-chromatin") }}

{{data_table("assets/base_asset_store/gwas/bmi/pulit_et_al_2018/analysis/pulit_et_al_2018_bmi_multi_tissue_chromatin_s_ldsc_cell_analysis_parquet_table.parquet",
id="multi-tissue-chromatin-bmi-table""
)}}


Again, by far the significant tissues are CNS-related, consistent with known biology.  There are also some miscellaneous significant tissue categories, like foreskin, mammary, and fetal stomach.  Whether these miscellaneous tissue results are biologically meaningful is unclear.