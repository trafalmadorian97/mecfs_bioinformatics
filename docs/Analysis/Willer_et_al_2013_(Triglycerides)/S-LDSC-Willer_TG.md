---
hide:
  - toc
---

# S-LDSC
I applied [Stratified Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md)[@finucane2018heritability] (S-LDSC) to summary statistics from Willer et al.'s GWAS of triglycerides[@global2013discovery] to identify possible key tissue and cell types affecting triglyceride levels.

## Reference Data Sources

I used the reference datasets recommended and preprocessed by the authors of the S-LSDC method[@finucane2018heritability].  These reference datasets are ultimately drawn from the following data sources:



- [The GTEx Project](../../Data_Sources/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../Data_Sources/Roadmap.md)
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../Data_Sources/Corces_et_al.md).
- The [ImmGen](../../Data_Sources/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset

## Results

### GTEx and Franke lab tissue expression data

I first applied S-LDSC to Willer et al.'s triglyceride GWAS using the Franke lab and GTEx datasets as a reference.  This produced  two hits, both in liver cells, which is consistent with the liver's known role as a hub of lipid metabolism. The results are shown below:


{{ plotly_embed("../../../_figs/willer_et_al_2023_tg_eur_multi_tissue_gene_expression_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="ldl-willer-sldsc-gtex", caption="Results of application of S-LDSC to Willer et al.'s triglyceride GWAS using the Franke lab/ GTEx dataset. Points are colored according to broad tissue category.  Large points correspond to cell/tissue types deemed significant by an application of the Benjamini-Hochberg procedure at an FDR of 0.01.") }}

### Roadmap Chromatin data

Next, I applied S-LDSC using the roadmap epigenetics dataset as a reference.  This produced no hits.


### ImmGen data

The next step is to use the S-LDSC reference data derived from the ImmGen project.  There were no significant cell types with this dataset.


### Corces ATAC-seq data

The next step is the use the Corces ATAC-seq data as a reference. Surprisingly, this produced a single strong hit in the monocyte cell type.

| Name    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------|--------------:|----------------------:|:--------------|
| Mono    |   4.402e-07   |           3.48457e-05 | True          |
| MEP     |   4.63353e-08 |           0.241071    | False         |
| GMP     |   2.94016e-08 |           0.312954    | False         |
| CMP     |   2.42517e-08 |           0.3478      | False         |
| Erythro |   4.58955e-08 |           0.351669    | False         |
| Bcell   |   1.96021e-08 |           0.379962    | False         |
| CLP     |   2.03106e-08 |           0.396499    | False         |
| MPP     |  -1.89262e-08 |           0.617103    | False         |
| HSC     |  -2.88261e-08 |           0.669921    | False         |
| NK      |  -2.96967e-08 |           0.69354     | False         |
| LMPP    |  -4.90885e-08 |           0.752546    | False         |
| CD4     |  -7.21404e-08 |           0.895904    | False         |
| CD8     |  -1.08379e-07 |           0.958275    | False         |⏎

The meaning of this hit is unclear.  It could be a proxy for a transcriptional program that is actually causal in another cell type, or it could reflect the fact that monocytes are macrophage precursors, and macrophages consume lipids to become foam cells when lipid levels are high[@steinberg2011cholesterol].


### Cahoy and GTEx-Brain data

The remaining two datasets pertain to the central nervous system.  There are no significant cell types with either of these datasets.  This is consistent with triglyceride levels being a largely non-neurological trait.
