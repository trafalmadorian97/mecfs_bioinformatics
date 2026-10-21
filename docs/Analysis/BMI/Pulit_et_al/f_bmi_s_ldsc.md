---
tags:
- S-LDSC
---
# S-LDSC

I applied [stratified linkage disequilibrium score regression](../../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) (S-LDSC)[@finucane2018heritability] to summary statistics from Pulit et al. meta-GWAS of body mass index (BMI)[@pulit2019meta].



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



The plot and table below show the results of S-LDSC applied to the BMI GWAS using the GTEx and Franke lab gene expression reference datasets. In the plot, the x-axis corresponds to cell type, while the y-axis corresponds to $-\log_{10}(p)$.  Points are colored according to tissue category.  Large points indicate cell/tissue types deemed significant by the Benjamini-Hochberg procedure at an FDR of 0.01[@benjamini1995controlling].

{{ plotly_embed("docs/_figs/pulit_et_al_2018_bmi_multi_tissue_gene_expression_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="bmi-sldsc-gene-expression") }}


{{
data_table("docs/_figs/pulit_et_al_2018_bmi_multi_tissue_gene_expression_s_ldsc_cell_analysis_parquet_table.parquet", id="bmi-gene-expression-table")
}}

Consistent with the known biology of obesity and our previous [MAGMA results](c_bmi_magma_gtex.md), all significant associations are in the central nervous system category.


### Roadmap/ ENCODE Chromatin data

I next applied S-LDSC using the Roadmap/ ENCODE epigenetic reference dataset.  The results are in the plot and expandable table below:

{{ plotly_embed("docs/_figs/pulit_et_al_2018_bmi_multi_tissue_chromatin_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="decode-bmi-sldsc-chromatin") }}

{{data_table("docs/_figs/pulit_et_al_2018_bmi_multi_tissue_chromatin_s_ldsc_cell_analysis_parquet_table.parquet",
id="multi-tissue-chromatin-bmi-table"
)}}


Again, by far the significant tissues are CNS-related, consistent with known biology.  There are also some miscellaneous significant tissue categories, like foreskin, mammary, and fetal stomach.  Whether these miscellaneous tissue results are biologically meaningful is unclear.

### ImmGen data


Next, I applied S-LDSC using the ImmGen project reference dataset.

The results are below

{{
data_table("docs/_figs/pulit_et_al_2018_bmi_immgen_s_ldsc_cell_analysis_parquet_table.parquet",id="bmi-immgen-table")
}}

There are no significant immune cell types.  This is consistent with BMI being primarily neural rather than immune.


### Corces et al. ATAC-seq data

The results of applying S-LDSC using the epigenetic reference data from the Corces et al. ATAC-seq analysis of blood cells are shown in the expandable table below.


{{
data_table("docs/_figs/pulit_et_al_2018_bmi_corces_atac_s_ldsc_cell_analysis_parquet_table.parquet",id="bmi-corces-table")
}}


Again, no-cell type are significant, consistent with BMI being not primarily driven by the immune system.




### Cahoy and GTEx-Brain data

The next two reference datasets pertain to the nervous system.


{{
data_table("docs/_figs/pulit_et_al_2018_bmi_cahoy_cns_s_ldsc_cell_analysis_parquet_table.parquet",id="bmi-cahoy-table")
}}


{{
data_table("docs/_figs/pulit_et_al_2018_bmi_gtex_brain_s_ldsc_cell_analysis_parquet_table.parquet",id="bmi-gtex-brain-table")
}}


There are no significant cell types from the Cahoy dataset.  There is a single significant cell type in the GTEx brain dataset: the cerebellum. Note that the S-LDSC annotations for brain cell types in the GTEx-Brain dataset are defined differently than the S-LDSC annotations for those same brain cell types in the [general GTEX dataset above](f_bmi_s_ldsc.md#gtex-and-franke-lab-tissue-expression-data): a genome region is assigned a cell-type annotation in the GTEx-Brain dataset only if it is preferentially associated with that brain region compared to other brain regions.  Thus the present results indicate BMI heritability enrichment in regions preferentially associated with the cerebellum.