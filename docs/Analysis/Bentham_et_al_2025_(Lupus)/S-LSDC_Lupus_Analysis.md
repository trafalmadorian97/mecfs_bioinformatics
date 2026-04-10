---
hide:
  - toc
---
# S-LDSC Analysis



I applied [Stratified Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md)[@finucane2018heritability] (S-LDSC) to summary statistics from Betham et al.'s GWAS of lupus[@bentham2015genetic].

## Reference Data Sources

I used the standard reference datasets prepared by the authors of the S-LDSC method[@finucane2018heritability].

- [The GTEx Project](../../Data_Sources/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../Data_Sources/Roadmap.md)
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../Data_Sources/Corces_et_al.md).
- The [ImmGen](../../Data_Sources/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset


## Results

### GTEx and Franke lab tissue expression data


The plot below shows the results of the application of S-LDSC to Bentham et al.'s study using the GTEx and Franke lab gene expression datasets. In the plot, the x-axis corresponds to cell type, while the y-axis corresponds to $-\log_{10}(p)$.  Points are colored according to broad tissue category.  Large points correspond to cell/tissue types deemed significant by an application of the Benjamini-Hochberg procedure at an FDR of 0.01[@benjamini1995controlling]. 



{{ plotly_embed("docs/_figs/bentham_2015_lupus_multi_tissue_gene_expression_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="lupus-sldsc-gtex", caption="") }}


As would be expected from an immunological disease like lupus, the top tissue types are all immune-related.


### Roadmap Chromatin data

I next applied S-LDSC using the reference dataset derived from the Roadmap epigenetic project.  The results are in the plot below:

{{ plotly_embed("docs/_figs/bentham_2015_lupus_multi_tissue_gene_expression_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="lupus-sldsc-chromatin", caption="") }}


Again, as would be expected of an immune condition, immune-related cell types rank highly.  Somewhat surprisingly, the signal for T-cells is stronger than the B-cell signal, despite lupus being considered a primarily immunoglobulin-driven disease.



### ImmGen data

Next, I applied S-LDSC using reference data from the ImmGen project.

There were no significant cell types.

The cell types with the lowest p values are given in the table below:


| Name                            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------------------------------|--------------:|----------------------:|:--------------|
| MF.11cloSer.Salm3.SI            |   1.11878e-07 |           0.0001879   | False         |
| DC.8-4-11b-.SLN                 |   1.19332e-07 |           0.000208158 | False         |
| NKT.44-NK1.1-.Th                |   9.40865e-08 |           0.000465922 | False         |
| DC.103+11b-.Lu                  |   9.82035e-08 |           0.000488146 | False         |
| DC.8-.Th                        |   1.11345e-07 |           0.000564316 | False         |
| DC.8+.Sp.ST                     |   8.85064e-08 |           0.00164973  | False         |
| DC.IIhilang-103-11b+.SLN        |   9.31931e-08 |           0.00178242  | False         |
| NKT.4-.Lv                       |   8.2538e-08  |           0.00211748  | False         |



Its surprising that no cell types were found to be significant in this analysis, given lupus's status as an immune-mediated disease.
To speculate, this may be due to differences between the mouse and human immune systems.

### Corces et al. ATAC-seq data

The results of applying S-LDSC using the epigenetic reference data from Corces et al. ATAC-seq analysis of blood cells are shown below.  


| Name    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------|--------------:|----------------------:|:--------------|
| Bcell   |   1.53276e-06 |           6.87243e-06 | True          |
| NK      |   1.38321e-06 |           0.00385197  | False         |
| CD8     |   1.26043e-06 |           0.0155016   | False         |
| CD4     |   1.02203e-06 |           0.023981    | False         |
| CLP     |   6.16263e-07 |           0.120099    | False         |
| Mono    |   4.66559e-07 |           0.247512    | False         |
| CMP     |   1.60724e-07 |           0.298654    | False         |
| MPP     |   1.39026e-07 |           0.331436    | False         |
| GMP     |   1.20103e-07 |           0.349802    | False         |
| LMPP    |   1.32611e-07 |           0.360171    | False         |
| MEP     |   9.50534e-08 |           0.373571    | False         |
| HSC     |   7.32517e-08 |           0.411087    | False         |
| Erythro |  -9.56772e-08 |           0.569965    | False         |⏎

There is one significant cell type: B-cells.  This is highly consistent with known status of lupus as a disease mediated by immunoglobulin.

### Cahoy and GTEx-Brain data

The next two reference datasets pertain to the nervous system.  The results of running S-LDSC with these two datasets are shown below:


| Name            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------|--------------:|----------------------:|:--------------|
| Oligodendrocyte |  -9.80353e-10 |              0.517435 | False         |
| Astrocyte       |  -8.08454e-09 |              0.637238 | False         |
| Neuron          |  -1.87689e-08 |              0.804622 | False         |⏎   


| Name                                    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------------------------------|--------------:|----------------------:|:--------------|
| Brain_Spinal_cord_(cervical_c-1)        |   5.77991e-08 |            0.00452981 | False         |
| Brain_Frontal_Cortex_(BA9)              |   1.90353e-08 |            0.139255   | False         |
| Brain_Substantia_nigra                  |   1.94392e-08 |            0.20111    | False         |
| Brain_Cerebellum                        |   9.35862e-09 |            0.344303   | False         |
| Brain_Cerebellar_Hemisphere             |   9.03005e-09 |            0.351129   | False         |
| Brain_Putamen_(basal_ganglia)           |   5.30034e-09 |            0.384532   | False         |
| Brain_Cortex                            |   4.89259e-09 |            0.405465   | False         |
| Brain_Amygdala                          |   3.51586e-09 |            0.42135    | False         |
| Brain_Hypothalamus                      |   1.05231e-09 |            0.478304   | False         |
| Brain_Anterior_cingulate_cortex_(BA24)  |  -2.08092e-09 |            0.548644   | False         |
| Brain_Nucleus_accumbens_(basal_ganglia) |  -1.57096e-08 |            0.791281   | False         |
| Brain_Hippocampus                       |  -1.64556e-08 |            0.804464   | False         |
| Brain_Caudate_(basal_ganglia)           |  -3.65828e-08 |            0.978487   | False         |


There are no significant cell types.  This is consistent with lupus being primarily non-neurological.



