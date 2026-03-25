---
hide:
  - toc
---
# S-LDSC Analysis



I applied [Stratified Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md)[@finucane2018heritability] (S-LDSC) to summary statistics from Nyeo et al.'s GWAS of blood EBV DNA levels[@nyeo2026population].

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

I first applied S-LDSC using the GTEx/Franke lab dataset as a reference.  At a false discovery rate of 0.01, there were no significant cell/tissue types.


### Roadmap epigenetic data


I next applied S-LDSC using the reference dataset derived from the Roadmap epigenetic project.  The results are in the plot below:

{{ plotly_embed("../../../_figs/nyeo_ebv_dna_multi_tissue_chromatin_cell_analysis_s_ldsc_plot/sldsc_scatter.html", id="ebv-sldsc-chromatin", caption="Results of application of S-LDSC to Nyeo et al.'s EBV DNA GWAS using the epigenetics reference dataset. Points are colored according to broad tissue category.  Large points correspond to cell/tissue types deemed significant by an application of the Benjamini-Hochberg procedure at an FDR of 0.01.") }}

At a false discovery rate of 0.01, there are a large number of significant cell types, most of which are T-cells or natural killer cells.  This is consistent with the idea that the key determinant of levels of EBV DNA is the ability of the immune system to contain a person's EBV infection, keeping it dormant.  

There is also one significant digestive system tissue type. This could be noise, or could perhaps reflect the presence of important immune cells in gut tissue.

That LDSC was able to pick up a much stronger signal from a chromatin epigenetics-based reference dataset compared to a gene-expression based reference dataset is interesting.  It is unclear why this would be the case.  Potentially, it means that EBV infections are controlled by genes whose differential expression levels are difficult to reliably measure.


### ImmGen data

Next, I applied S-LDSC using reference data from the ImmGen project.

There were no significant cell types.

The cell types with the lowest p values are given in the table below:



| Name                            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------------------------------|--------------:|----------------------:|:--------------|
| T.8SP24-.Th                     |   3.082e-09   |           0.000625936 | False         |
| NKT.4-.Lv                       |   2.58537e-09 |           0.000767295 | False         |
| T.4.PLN.BDC                     |   3.20531e-09 |           0.000847326 | False         |
| T.8Nve.MLN                      |   2.45397e-09 |           0.00110938  | False         |
| NKT.4+.Lv                       |   2.11435e-09 |           0.00114223  | False         |
| T.8Nve.LN                       |   2.6488e-09  |           0.00151794  | False         |
| T.4FP3+25+.Sp                   |   2.31043e-09 |           0.0015654   | False         |
| T.4SP69+.Th                     |   2.52046e-09 |           0.00255811  | False         |
| T.4SP24-.Th                     |   2.57318e-09 |           0.00291993  | False         |
| Tgd.vg2-.act.Sp                 |   1.95198e-09 |           0.0031735   | False         |
| T.4Mem.Sp                       |   3.69048e-09 |           0.0036837   | False         |
| Tgd.vg2+.act.Sp                 |   1.70996e-09 |           0.00368639  | False         |
| T.4FP3-.Sp                      |   2.56193e-09 |           0.00436498  | False         |
| T.4.LN.BDC                      |   2.51854e-09 |           0.00500759  | False         |
| T.4.Pa.BDC                      |   3.65641e-09 |           0.00502588  | False         |
| T.8SP69+.Th                     |   2.60463e-09 |           0.00506605  | False         |
| NKT.4-.Sp                       |   1.64369e-09 |           0.00586899  | False         |
| NKT.44+NK1.1+.Th                |   2.07751e-09 |           0.00631805  | False         |
| NKT.44-NK1.1-.Th                |   2.1653e-09  |           0.0063773   | False         |
| T.4SP24int.Th                   |   2.19929e-09 |           0.00648667  | False         |
| Tgd.vg2-.Sp                     |   2.03306e-09 |           0.00652517  | False         |
| Tgd.vg2+.Sp                     |   1.86698e-09 |           0.00765776  | False         |
| T.8Mem.Sp                       |   3.51069e-09 |           0.00781653  | False         |
| T.8Mem.Sp.OT1.d100.LisOva       |   3.96871e-09 |           0.00816987  | False         |
| Tgd.vg5-.act.IEL                |   3.31278e-09 |           0.00946079  | False         |


Consistent with the results above, the lowest p values are all in T or natural killer cells.

It is unclear why the signal from this reference dataset is not sufficient to produce significant cell types.

This could be related to:

- Differences between the mouse and human immune system
- Problems with measuring gene expression levels of key genes
- Since ImmGen is an immune dataset, for a cell type to be significant, there must a differential signal in that cell type compared to other immune cell types.  If the EBV DNA GWAS instead provides a generalized immune signal, this could result in no significant cell types.


### Corces et al. ATAC-seq data

The results of applying S-LDSC using the epigenetic reference data from Corces et al. ATAC-seq analysis of blood cells are shown below.  


| Name    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------|--------------:|----------------------:|:--------------|
| CD8     |   7.49151e-08 |           6.86179e-07 | True          |
| CD4     |   6.03118e-08 |           4.32505e-05 | True          |
| NK      |   6.50971e-08 |           0.000225251 | True          |
| Bcell   |   3.61798e-08 |           0.000610793 | True          |
| Mono    |   2.97371e-08 |           0.0150854   | False         |
| GMP     |   1.45912e-08 |           0.0519185   | False         |
| HSC     |   1.33929e-08 |           0.0645347   | False         |
| LMPP    |   1.44259e-08 |           0.0895894   | False         |
| MPP     |   1.05311e-08 |           0.0980141   | False         |
| CLP     |   1.79604e-08 |           0.100229    | False         |
| CMP     |   7.3202e-09  |           0.17619     | False         |
| MEP     |   3.95883e-09 |           0.308036    | False         |
| Erythro |   7.45506e-10 |           0.475241    | False         |


T-cells (CD4 and CD8), natural killer (NK) cells, and B-cells are all significant.   Again, this is consistent with an immune mechanism of control of EBV levels.  Also, this is yet another line of evidence reinforcing the primacy of T-cells as the key EBV-level-determining cells.



### Cahoy and GTEx-Brain data

The next two reference datasets pertain to the nervous system.  The results of running S-LDSC with these two datasets are shown below:

| Name            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------|--------------:|----------------------:|:--------------|
| Astrocyte       |   1.32785e-09 |             0.0924435 | False         |
| Neuron          |  -6.41258e-10 |             0.828212  | False         |
| Oligodendrocyte |  -1.82847e-09 |             0.980079  | False         |


| Name                                    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------------------------------|--------------:|----------------------:|:--------------|
| Brain_Spinal_cord_(cervical_c-1)        |   1.15863e-09 |             0.0776225 | False         |
| Brain_Amygdala                          |   9.11299e-10 |             0.139974  | False         |
| Brain_Frontal_Cortex_(BA9)              |   9.32179e-10 |             0.144977  | False         |
| Brain_Cerebellar_Hemisphere             |   8.02938e-10 |             0.15864   | False         |
| Brain_Cerebellum                        |   5.50179e-10 |             0.237984  | False         |
| Brain_Hypothalamus                      |  -2.07267e-10 |             0.690906  | False         |
| Brain_Anterior_cingulate_cortex_(BA24)  |  -2.61295e-10 |             0.704288  | False         |
| Brain_Nucleus_accumbens_(basal_ganglia) |  -3.27695e-10 |             0.741726  | False         |
| Brain_Putamen_(basal_ganglia)           |  -4.7383e-10  |             0.748924  | False         |
| Brain_Substantia_nigra                  |  -5.81302e-10 |             0.859182  | False         |
| Brain_Hippocampus                       |  -6.91034e-10 |             0.867084  | False         |
| Brain_Cortex                            |  -6.85036e-10 |             0.898975  | False         |
| Brain_Caudate_(basal_ganglia)           |  -8.09745e-10 |             0.914915  | False         |


There are no significant cell types, which is consistent with the determinants of EBV DNA levels being primarily non-neurological.


