---
hide:
- navigation
- toc
---
# S-LDSC Analysis of SCH



I used [Stratified Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) (S-LDSC) to analyze summary statistics from a GWAS of Schizophrenia from the Psychiatric Genetic Consortium[@trubetskoy2022mapping].  For the current analysis, I focused only European participants.


I used the standard reference data sources to associate chromatin regions with cell or tissue types. There reference data sources were:

- [The GTEx Project](../../Data_Sources/GTEx_Project/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../Data_Sources/Roadmap_Epigenetic_Project/Roadmap.md)
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../Data_Sources/Corces_Et_Al_Blood_Cell_Epigenetics/Corces_et_al.md).
- The [ImmGen](../../Data_Sources/Immgen_Project/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset

## Results

### GTEx and Franke lab tissue expression data
The plot below illustrates the pattern of S-LDSC p-values across cell types. The table shows the cell types with the lowest p-values.


![s-ldsc-sch-gene-exp](https://github.com/user-attachments/assets/885f815c-2934-486e-aef8-0d3c799a1c5d)


| Name                                                           |   Coefficient |   Coefficient_P_value | Reject Null   |
|:---------------------------------------------------------------|--------------:|----------------------:|:--------------|
| A08.186.211.730.885.287.500.Cerebral.Cortex                    |   8.47046e-08 |           3.98912e-19 | True          |
| A08.186.211.464.710.225.Entorhinal.Cortex                      |   8.82079e-08 |           5.25314e-19 | True          |
| A08.186.211.464.Limbic.System                                  |   8.59402e-08 |           5.35036e-19 | True          |
| Brain_Frontal_Cortex_(BA9)                                     |   7.50537e-08 |           2.70948e-18 | True          |
| Brain_Cortex                                                   |   7.43564e-08 |           3.53839e-18 | True          |
| Brain_Anterior_cingulate_cortex_(BA24)                         |   7.3208e-08  |           1.8928e-16  | True          |
| A08.186.211.464.405.Hippocampus                                |   7.93705e-08 |           2.02728e-16 | True          |
| A08.186.211.Brain                                              |   8.6903e-08  |           2.41997e-16 | True          |
| Brain_Hippocampus                                              |   6.32123e-08 |           8.39121e-15 | True          |
| Brain_Cerebellar_Hemisphere                                    |   5.94796e-08 |           7.06571e-14 | True          |
| Brain_Amygdala                                                 |   6.0501e-08  |           3.274e-13   | True          |
| Brain_Putamen_(basal_ganglia)                                  |   6.33216e-08 |           1.64665e-12 | True          |
| A08.186.211.132.Brain.Stem                                     |   6.74624e-08 |           2.27026e-12 | True          |
| Brain_Cerebellum                                               |   5.66552e-08 |           4.05769e-12 | True          |
| A08.186.211.730.885.287.500.571.735.Visual.Cortex              |   6.7991e-08  |           8.32935e-12 | True          |
| Brain_Nucleus_accumbens_(basal_ganglia)                        |   5.99918e-08 |           1.19195e-11 | True          |
| Brain_Caudate_(basal_ganglia)                                  |   6.11288e-08 |           2.41973e-11 | True          |
| A08.186.211.730.885.287.500.270.Frontal.Lobe                   |   6.86582e-08 |           1.13253e-10 | True          |
| Brain_Substantia_nigra                                         |   5.53578e-08 |           1.26689e-10 | True          |
| A08.186.211.730.885.287.500.670.Parietal.Lobe                  |   6.54824e-08 |           1.3482e-10  | True          |
| A08.186.211.132.810.428.200.Cerebellum                         |   6.13038e-08 |           4.10792e-10 | True          |
| Brain_Hypothalamus                                             |   4.90558e-08 |           1.61256e-09 | True          |
| A08.186.211.865.428.Metencephalon                              |   5.9104e-08  |           3.64513e-09 | True          |
| A08.186.211.730.885.287.249.Basal.Ganglia                      |   5.14062e-08 |           9.02497e-07 | True          |
| Brain_Spinal_cord_(cervical_c-1)                               |   4.00887e-08 |           2.68969e-06 | True          |
| A08.186.211.730.317.Diencephalon                               |   5.77876e-08 |           5.75561e-06 | True          |
| A08.186.211.653.Mesencephalon                                  |   3.83026e-08 |           8.03001e-05 | True          |
| A08.186.211.730.317.357.Hypothalamus                           |   5.184e-08   |           0.000119281 | True          |
| A08.186.211.730.885.287.249.487.Corpus.Striatum                |   3.79874e-08 |           0.000125097 | True          |
| A09.371.729.Retina                                             |   3.16838e-08 |           0.000147619 | True          |
| A11.872.653.Neural.Stem.Cells                                  |   4.28784e-08 |           0.00020438  | True          |
| A08.186.211.730.317.357.352.435.Hypothalamo.Hypophyseal.System |   4.51566e-08 |           0.00146605  | True          |
| A11.118.637.555.567.562.440.Precursor.Cells..B.Lymphoid        |   3.17402e-08 |           0.00412467  | False         |
| A11.118.637.555.567.569.T.Lymphocytes                          |   3.0355e-08  |           0.00699123  | False         |
| A15.145.229.637.555.567.569.200.CD4.Positive.T.Lymphocytes     |   2.72466e-08 |           0.0118845   | False         |
| A15.145.Blood                                                  |   2.7143e-08  |           0.022607    | False         |
| A11.118.637.555.567.569.200.700.T.Lymphocytes..Regulatory      |   1.71779e-08 |           0.060104    | False         |
| A15.145.229.Blood.Cells                                        |   2.08363e-08 |           0.0640231   | False         |
| A11.872.378.590.635.Granulocyte.Macrophage.Progenitor.Cells    |   1.63959e-08 |           0.0854463   | False         |
| A15.382.490.315.583.Neutrophils                                |   1.63362e-08 |           0.10354     | False         |
| A11.118.637.Leukocytes                                         |   1.58603e-08 |           0.115801    | False         |
| Muscle_Skeletal                                                |   1.10094e-08 |           0.129756    | False         |
| Artery_Tibial                                                  |   1.03624e-08 |           0.176228    | False         |
| A15.378.316.580.Monocytes                                      |   1.13324e-08 |           0.199719    | False         |
| A15.382.520.604.800.Palatine.Tonsil                            |   8.02661e-09 |           0.223744    | False         |
| Pituitary                                                      |   7.09361e-09 |           0.236527    | False         |
| A15.145.229.637.555.Leukocytes..Mononuclear                    |   5.93882e-09 |           0.28341     | False         |
| A15.382.490.555.567.Lymphocytes                                |   6.12958e-09 |           0.299802    | False         |


The specificity of the signal for central nervous system cell types is consistent with the known biology of schizophrenia as neurological condition.  Of the non-significant cell types, immune cells are closest to being significant, which is consistent with theories of the role of the immune system in the development of Schizophrenia  (See pg. 50 and pg. 1647 in Kandel et al.[@kandel2021principles]).



### Roadmap Chromatin data

I next applied S-LDSC to the schizophrenia GWAS using reference data generated by Finucane et al. from the [Roadmap Epigenetics Project](../../Data_Sources/Roadmap_Epigenetic_Project/Roadmap.md)

The results are shown in the graph and table below

![s-ldsc-sch-chromatin](https://github.com/user-attachments/assets/a0449b42-9bc3-4347-82d3-bee9eec1f88c)

| Name                                                                     |   Coefficient |   Coefficient_P_value | Reject Null   |
|:-------------------------------------------------------------------------|--------------:|----------------------:|:--------------|
| Fetal_Brain_Male__DNase                                                  |   1.68738e-06 |           4.74122e-22 | True          |
| Fetal_Brain_Female__DNase                                                |   1.72679e-06 |           2.12598e-19 | True          |
| Fetal_Brain_Female__H3K4me3                                              |   2.77978e-06 |           4.96519e-17 | True          |
| Fetal_Brain_Male__H3K4me1                                                |   7.68472e-07 |           1.65669e-14 | True          |
| Brain_Germinal_Matrix__H3K4me3                                           |   2.7504e-06  |           2.05815e-14 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K27ac                            |   9.67876e-07 |           6.35169e-14 | True          |
| Brain_Anterior_Caudate__H3K4me1                                          |   6.73639e-07 |           2.08363e-12 | True          |
| Brain_Inferior_Temporal_Lobe__H3K27ac                                    |   6.39341e-07 |           2.33773e-12 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K4me3                            |   2.8196e-06  |           2.69558e-12 | True          |
| Brain_Angular_Gyrus__H3K4me3                                             |   2.71607e-06 |           5.17278e-12 | True          |
| Brain_Angular_Gyrus__H3K4me1                                             |   8.00012e-07 |           9.72589e-12 | True          |
| Brain_Inferior_Temporal_Lobe__H3K4me3                                    |   2.19647e-06 |           1.18606e-11 | True          |
| Brain_Anterior_Caudate__H3K27ac                                          |   5.85791e-07 |           1.32123e-11 | True          |
| Brain_Angular_Gyrus__H3K27ac                                             |   8.13771e-07 |           1.40954e-11 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K4me1                            |   1.01035e-06 |           4.05928e-11 | True          |
| Fetal_Brain_Female__H3K4me1                                              |   1.05581e-06 |           1.26443e-10 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K9ac                             |   2.26479e-06 |           3.4872e-10  | True          |
| Brain_Anterior_Caudate__H3K9ac                                           |   1.0844e-06  |           4.82425e-10 | True          |
| Fetal_Brain_Female__H3K36me3                                             |   5.79725e-07 |           7.90235e-10 | True          |
| Brain_Cingulate_Gyrus__H3K9ac                                            |   1.47375e-06 |           1.69795e-09 | True          |
| Brain_Angular_Gyrus__H3K9ac                                              |   1.38489e-06 |           5.17934e-09 | True          |
| Brain_Cingulate_Gyrus__H3K4me1                                           |   5.86637e-07 |           5.19947e-09 | True          |
| Ganglion_Eminence_derived_primary_cultured_neurospheres__H3K4me3         |   1.57836e-06 |           5.34906e-09 | True          |
| Brain_Cingulate_Gyrus__H3K4me3                                           |   2.03552e-06 |           8.42614e-09 | True          |
| Brain_Anterior_Caudate__H3K4me3                                          |   1.5262e-06  |           2.21405e-08 | True          |
| Brain_Germinal_Matrix__H3K4me1                                           |   9.61141e-07 |           2.70825e-08 | True          |
| Brain_Germinal_Matrix__H3K36me3                                          |   7.31311e-07 |           3.62691e-08 | True          |
| Cortex_derived_primary_cultured_neurospheres__H3K4me3                    |   1.9088e-06  |           7.93075e-08 | True          |
| Brain_Cingulate_Gyrus__H3K27ac                                           |   5.22876e-07 |           1.20097e-07 | True          |
| Brain_Hippocampus_Middle__H3K4me1                                        |   3.88849e-07 |           1.44024e-07 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K36me3                           |   7.34169e-07 |           8.23128e-07 | True          |
| Brain_Inferior_Temporal_Lobe__H3K9ac                                     |   9.45094e-07 |           1.01577e-06 | True          |
| Brain_Angular_Gyrus__H3K36me3                                            |   7.09494e-07 |           1.26005e-06 | True          |
| Brain_Inferior_Temporal_Lobe__H3K4me1                                    |   6.24768e-07 |           2.82109e-06 | True          |
| Brain_Hippocampus_Middle__H3K4me3                                        |   9.85198e-07 |           1.80801e-05 | True          |
| Ganglion_Eminence_derived_primary_cultured_neurospheres__H3K4me1         |   5.15632e-07 |           2.35424e-05 | True          |
| Brain_Inferior_Temporal_Lobe__H3K36me3                                   |   4.74221e-07 |           2.69696e-05 | True          |
| Brain_Hippocampus_Middle__H3K27ac                                        |   3.29051e-07 |           6.20309e-05 | True          |
| Primary_T_cells_from_cord_blood__DNase                                   |   7.41953e-07 |           0.000137651 | True          |
| Primary_B_cells_from_peripheral_blood__DNase                             |   6.10191e-07 |           0.000146033 | True          |
| Brain_Substantia_Nigra__H3K4me1                                          |   3.04377e-07 |           0.000193331 | True          |
| Brain_Anterior_Caudate__H3K36me3                                         |   4.48397e-07 |           0.000217432 | True          |
| Cortex_derived_primary_cultured_neurospheres__H3K36me3                   |   4.0351e-07  |           0.000398369 | True          |
| Spleen_ENTEX__H3K27ac                                                    |   1.54866e-07 |           0.000520367 | True          |
| Ganglion_Eminence_derived_primary_cultured_neurospheres__H3K36me3        |   3.65732e-07 |           0.000644606 | True          |
| Brain_Hippocampus_Middle__H3K36me3                                       |   3.00257e-07 |           0.000730833 | True          |
| Cortex_derived_primary_cultured_neurospheres__H3K4me1                    |   2.95148e-07 |           0.000907247 | True          |
| Brain_Cingulate_Gyrus__H3K36me3                                          |   4.24078e-07 |           0.00107608  | False         |
| Pancreatic_Islets__H3K4me1                                               |   2.63705e-07 |           0.00111471  | False         |
| Brain_Substantia_Nigra__H3K9ac                                           |   5.91492e-07 |           0.00140614  | False         |
| Esoph-GJ_ENTEX__H3K4me3                                                  |   7.76109e-07 |           0.0014798   | False         |
| Heart-Atrial_ENTEX__H3K27ac                                              |   1.7608e-07  |           0.0020061   | False         |
| skeletal_muscle_ENTEX__H3K4me1                                           |   1.2165e-07  |           0.00232854  | False         |


As we saw with the gene expression data, signal is highly concentrated in central nervous system tissues.  Unlike with the gene expression data, with the roadmap data we also seem some immune-related tissues and cell types reaching significance, consistent with the purposed role of the immune system in schizophrenia pathogenesis[@kandel2021principles].


### Immgen Data

The next step is to use the SLDSC reference data derived from the Immgen project.
In this case, no cell types are significant. The top non-significant cell types are shown below:

| Name                            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------------------------------|--------------:|----------------------:|:--------------|
| B.FrE.BM                        |   3.69752e-08 |           0.000926076 | False         |
| preT.ETP-2A.Th                  |   3.41931e-08 |           0.0023043   | False         |
| NKT.44-NK1.1-.Th                |   3.05201e-08 |           0.00344786  | False         |
| B.T1.Sp                         |   2.71268e-08 |           0.00707957  | False         |
| CD19Control                     |   2.72862e-08 |           0.00819193  | False         |
| T.4int8+.Th                     |   2.81369e-08 |           0.00889715  | False         |
| B.Fo.Sp                         |   2.81622e-08 |           0.00942477  | False         |
| proB.FrA.BM                     |   2.62878e-08 |           0.01131     | False         |


This inability to localize the signal to a specific immune cell type using the immgen data may reflect insufficient power in the schhizophrenia GWAS, or a weakness of the mouse-derived immgen data as a reference for S-LDSC.


### Corces et al. ATAC-seq data

Similarly, we are not able to identify any specific cell type using the Corces ATAC-seq dataset:

| Name    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------|--------------:|----------------------:|:--------------|
| Bcell   |   3.05376e-07 |             0.0179416 | False         |
| CD8     |   1.30937e-07 |             0.180235  | False         |
| Erythro |   4.99447e-08 |             0.399927  | False         |
| NK      |   2.49454e-08 |             0.426503  | False         |
| CD4     |   3.11632e-09 |             0.490951  | False         |
| MEP     |  -1.24416e-08 |             0.542468  | False         |
| Mono    |  -6.1538e-08  |             0.610702  | False         |
| CMP     |  -9.74273e-08 |             0.797071  | False         |


### Cahoy and GTEx-Brain data

A expected, the Cahoy dataset strongly localizes the Schizophrenia GWAS signal in neurons:


| Name    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------|--------------:|----------------------:|:--------------|
| Neuron          |   4.48721e-08 |           5.83951e-06 | True          |
| Oligodendrocyte |   1.91268e-08 |           0.0638051   | False         |
| Astrocyte       |  -1.22911e-08 |           0.872973    | False         |⏎


The GTEx brain dataset points to the cortex as a key site of disease activity. Here are tissue types with the lowest p-values:



| Name                                    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------------------------------|--------------:|----------------------:|:--------------|
| Brain_Frontal_Cortex_(BA9)              |   4.66646e-08 |           1.63948e-09 | True          |
| Brain_Cortex                            |   5.20674e-08 |           2.42723e-09 | True          |
| Brain_Anterior_cingulate_cortex_(BA24)  |   3.91227e-08 |           3.47957e-06 | True          |
| Brain_Cerebellar_Hemisphere             |   2.96201e-08 |           0.00209779  | True          |
| Brain_Cerebellum                        |   1.97362e-08 |           0.0162309   | False         |
| Brain_Putamen_(basal_ganglia)           |   5.0827e-10  |           0.479543    | False         |
| Brain_Nucleus_accumbens_(basal_ganglia) |  -5.24356e-10 |           0.523653    | False         |
| Brain_Hippocampus                       |  -2.20603e-09 |           0.584205    | False         |
| Brain_Amygdala                          |  -2.75695e-09 |           0.623402    | False         |
| Brain_Caudate_(basal_ganglia)           |  -4.64408e-09 |           0.697234    | False         |
| Brain_Hypothalamus                      |  -2.18885e-08 |           0.993308    | False         |
| Brain_Spinal_cord_(cervical_c-1)        |  -2.66567e-08 |           0.996263    | False         |
| Brain_Substantia_nigra                  |  -3.35068e-08 |           0.999814    | False         |⏎


## Reproducing

Reproduce these results using the [Schizophnreia Analysis Script][mecfs_bio.analysis.schizophrenia_analysis].