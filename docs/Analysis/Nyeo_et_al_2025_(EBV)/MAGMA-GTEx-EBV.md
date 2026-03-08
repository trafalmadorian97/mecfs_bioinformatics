---
hide:
  - toc
---
# MAGMA GTEx


I applied [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to EBV DNA GWAS of Nyeo et al.[@nyeo2026population] using bulk RNAseq data from [GTEx](../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex] as a reference. The gene-property analysis results are plotted below.


{{ plotly_embed("../../../_figs/nyeo_ebv_dna_magma_bar_plot/magma_gene_set_plot.html", id="ebv-gtex-magma", caption="Result of MAGMA GTEx applied to Nyeo et al.'s EBV GWAS. The x axis corresponds to tissue type, while the y axis measures MAGMA significance level.") }}

The four significant tissue/cell types are _Ileum Lymphode Aggregate, Whole Blood, Spleen, and EBV-Transformed Lymphocytes_. All are immune-related. Given that the primary determinant of a person's EBV DNA levels is the ability of their immune system to keep their EBV infection under control and in a dormant state, this makes sense.