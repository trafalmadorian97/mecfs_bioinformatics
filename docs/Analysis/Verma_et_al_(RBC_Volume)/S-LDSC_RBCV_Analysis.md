---
hide:
    - navigation
    - toc
---

# S-LDSC RBCV Analysis

I applied [Stratified Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) (S-LDSC) to summary statistics from a GWAS of Red Blood Cell Volume from the [Million Veterans Program](../../Data_Sources/Million_Veterans_Program.md).

I used standard S-LDSC reference datasets. These reference datasets were originally derived from:

- [The GTEx Project](../../Data_Sources/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../Data_Sources/Roadmap.md)
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../Data_Sources/Corces_et_al.md).
- The [ImmGen](../../Data_Sources/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset


## Results

### GTEx and Franke lab tissue expression data
The plot below illustrates the pattern of S-LDSC p-values across cell types. The table shows the cell types with the lowest p-values.


![s-ldsc_rcvol_expression](https://github.com/user-attachments/assets/2bef69a9-7f86-45fc-b42c-894bb00cfafe)

| Name                                                           |   Coefficient |   Coefficient_P_value | Reject Null   |
|:---------------------------------------------------------------|--------------:|----------------------:|:--------------|
| A15.145.300.Fetal.Blood                                        |   8.02657e-08 |           6.59384e-09 | True          |
| A15.145.Blood                                                  |   5.94179e-08 |           7.3992e-06  | True          |
| A11.443.Erythroid.Cells                                        |   4.41781e-08 |           1.19728e-05 | True          |
| A11.627.635.Myeloid.Progenitor.Cells                           |   4.62524e-08 |           0.000107559 | True          |
| A11.872.378.Hematopoietic.Stem.Cells                           |   4.41011e-08 |           0.000172671 | True          |
| A15.382.Immune.System                                          |   3.99963e-08 |           0.000342322 | False         |
| Whole_Blood                                                    |   3.64886e-08 |           0.000399347 | False         |
| A15.145.229.188.Blood.Platelets                                |   3.14462e-08 |           0.000464739 | False         |
| A11.872.378.590.817.Megakaryocyte.Erythroid.Progenitor.Cells   |   3.67321e-08 |           0.000603583 | False         |
| A15.145.229.Blood.Cells                                        |   4.13485e-08 |           0.000654442 | False         |
| A11.118.637.555.567.562.B.Lymphocytes                          |   3.57229e-08 |           0.000800263 | False         |
| Spleen                                                         |   3.89864e-08 |           0.000854701 | False         |
| A11.118.637.555.567.569.T.Lymphocytes                          |   2.90771e-08 |           0.00104616  | False         |
| Colon_Transverse                                               |   2.6992e-08  |           0.00126636  | False         |
| A15.145.229.637.555.567.569.200.CD4.Positive.T.Lymphocytes     |   2.53835e-08 |           0.00230767  | False         |
| A15.382.490.555.567.537.Killer.Cells..Natural                  |   3.36356e-08 |           0.00263635  | False         |
| A11.627.340.360.Granulocyte.Precursor.Cells                    |   2.73892e-08 |           0.00301885  | False         |
| A15.145.229.637.555.567.562.725.Plasma.Cells                   |   3.18685e-08 |           0.0035067   | False         |
| A11.118.637.Leukocytes                                         |   3.64806e-08 |           0.00398084  | False         |
| A15.382.520.604.800.Palatine.Tonsil                            |   2.93966e-08 |           0.00419334  | False         |
| A11.118.637.555.567.569.200.700.T.Lymphocytes..Regulatory      |   2.38386e-08 |           0.00485846  | False         |
| A11.872.Stem.Cells                                             |   2.63215e-08 |           0.00487292  | False         |
| A11.118.637.555.567.562.440.Precursor.Cells..B.Lymphoid        |   2.5128e-08  |           0.00587332  | False         |
| A11.872.378.590.635.Granulocyte.Macrophage.Progenitor.Cells    |   2.81929e-08 |           0.00629756  | False         |
| A03.556.124.684.Intestine..Small                               |   2.71434e-08 |           0.0065321   | False         |
| A15.382.490.555.567.Lymphocytes                                |   3.1885e-08  |           0.00745577  | False         |
| A03.556.249.124.Ileum                                          |   2.81618e-08 |           0.00763379  | False         |
| Cells_EBV-transformed_lymphocytes                              |   2.55655e-08 |           0.00789262  | False         |

The 5 significant cell types are: "A15.145.300.Fetal.Blood", "A15.145.Blood", "A11.443.Erythroid.Cells", "A11.627.635.Myeloid.Progenitor.Cells".

These findings makes sense, and are consistent with known biology: "Erethryoid cells" are a class that includes standard red blood cells (erythrocytes), while myeloid progenitor cells are ancestors to red blood cells (see the figure below). 


![haemopoetic-cells](https://github.com/user-attachments/assets/ccb67146-5c88-4168-be81-769a55e6b844)

Figure: By A. Rad and M. Häggström. CC-BY-SA 3.0 license

### Roadmap Chromatin data

I next applied S-LDSC to the red blood cell volume GWAS using reference data generated by Finucane et al. from the [Roadmap Epigenetics Project](../../Data_Sources/Roadmap.md).

The following graph and table show the results:

![s-lsdc-rbcvol-chromatin](https://github.com/user-attachments/assets/8b21bd47-daf3-4d52-bc51-70f8a3d10946)

| Name                                                                     |   Coefficient |   Coefficient_P_value | Reject Null   |
|:-------------------------------------------------------------------------|--------------:|----------------------:|:--------------|
| Primary_hematopoietic_stem_cells_G-CSF-mobilized_Female__H3K27ac         |   9.7986e-07  |           4.30062e-07 | True          |
| Primary_hematopoietic_stem_cells_G-CSF-mobilized_Female__H3K4me1         |   3.89211e-07 |           7.36713e-07 | True          |
| Primary_hematopoietic_stem_cells_G-CSF-mobilized_Female__DNase           |   9.67532e-07 |           1.22959e-06 | True          |
| Primary_hematopoietic_stem_cells__H3K4me1                                |   7.08867e-07 |           3.25765e-06 | True          |
| Primary_hematopoietic_stem_cells_G-CSF-mobilized_Male__H3K4me1           |   4.06036e-07 |           3.37892e-06 | True          |
| Primary_T_killer_memory_cells_from_peripheral_blood__H3K27ac             |   8.61697e-07 |           3.55327e-06 | True          |
| Primary_hematopoietic_stem_cells_short_term_culture__H3K4me3             |   1.833e-06   |           5.35886e-06 | True          |
| Primary_hematopoietic_stem_cells_G-CSF-mobilized_Male__DNase             |   8.90935e-07 |           7.23909e-06 | True          |
| Primary_T_helper_memory_cells_from_peripheral_blood_2__H3K27ac           |   6.6463e-07  |           2.52533e-05 | True          |
| Primary_hematopoietic_stem_cells_G-CSF-mobilized_Male__H3K4me3           |   1.0486e-06  |           2.63718e-05 | True          |
| Spleen_ENTEX__H3K27ac                                                    |   2.05768e-07 |           4.87077e-05 | True          |
| Fetal_Thymus__H3K36me3                                                   |   2.79975e-07 |           5.17104e-05 | True          |
| Primary_hematopoietic_stem_cells_short_term_culture__H3K4me1             |   3.5475e-07  |           8.47972e-05 | True          |
| Primary_hematopoietic_stem_cells_G-CSF-mobilized_Female__H3K4me3         |   1.75499e-06 |           0.000127767 | True          |
| Primary_hematopoietic_stem_cells_short_term_culture__H3K36me3            |   3.44274e-07 |           0.000129613 | True          |
| Spleen_ENTEX__H3K4me3                                                    |   3.49278e-07 |           0.000195307 | True          |
| Primary_neutrophils_from_peripheral_blood__H3K4me3                       |   6.33005e-07 |           0.000318123 | True          |
| Primary_B_cells_from_peripheral_blood__DNase                             |   4.68288e-07 |           0.000512658 | False         |
| Primary_neutrophils_from_peripheral_blood__H3K4me1                       |   1.87677e-07 |           0.000658363 | False         |
| Spleen_ENTEX__H3K4me1                                                    |   9.59291e-07 |           0.00107238  | False         |
| Rectal_Mucosa_Donor_31__H3K9ac                                           |   5.6792e-07  |           0.00116076  | False         |
| Stomach_ENTEX__H3K4me3                                                   |   7.71536e-07 |           0.00122671  | False         |
| SI-Term-Ileum_ENTEX__H3K4me3                                             |   6.86199e-07 |           0.00140399  | False         |
| Primary_hematopoietic_stem_cells__H3K4me3                                |   1.4126e-06  |           0.00167475  | False         |
| Primary_T_helper_cells_from_peripheral_blood__H3K4me3                    |   7.52953e-07 |           0.00177749  | False         |
| Fetal_Thymus__DNase                                                      |   4.10651e-07 |           0.00233247  | False         |
| SI-Term-Ileum_ENTEX__H3K27ac                                             |   2.38673e-07 |           0.0023326   | False         |
| Primary_Natural_Killer_cells_from_peripheral_blood__H3K4me1              |   2.12516e-07 |           0.00384186  | False         |


Again, consistent with known biology, the significant cell types are all in the blood/immune category and related to or descended from haematopoietic stem cells.  

The presence of lymphocytes like T-cells on the list of significant cell types is interesting.  My guess is that this may reflect epigenetic features shared by all cells in haematopoietic stem cell lineage, rather than true etiological role for lymphocytes in red blood cell volume.

### Imgen Data

The next step is to use the SLDSC reference data derived from the Immgen project.

![immgen-rbc-plot](https://github.com/user-attachments/assets/a227056b-2b3c-4a61-99a3-3717a6e23dca)

| Name    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------------------------------|--------------:|----------------------:|:--------------|
| preT.DN2.Th                     |   7.01741e-08 |           7.55644e-08 | True          |
| T.8SP69+.Th                     |   6.10773e-08 |           4.08943e-07 | True          |
| preT.ETP-2A.Th                  |   6.13855e-08 |           4.64764e-07 | True          |
| preT.DN2-3.Th                   |   5.88185e-08 |           5.87392e-07 | True          |
| T.4int8+.Th                     |   5.54558e-08 |           6.74419e-07 | True          |
| preT.DN2A.Th                    |   5.66926e-08 |           1.21757e-06 | True          |
| Tgd.vg2+24ahi.Th                |   5.37253e-08 |           2.51879e-06 | True          |
| CD4.96h.LN                      |   5.33706e-08 |           4.68263e-06 | True          |
| T.8SP24int.Th                   |   4.98821e-08 |           4.85374e-06 | True          |
| Tgd.Th                          |   5.3534e-08  |           6.5081e-06  | True          |
| preT.DN3-4.Th                   |   5.45315e-08 |           7.09887e-06 | True          |
| preT.DN3B.Th                    |   5.52977e-08 |           1.07887e-05 | True          |
| SC.LT34F.BM                     |   4.19728e-08 |           1.54079e-05 | True          |
| Tgd.vg2-.Sp                     |   4.27974e-08 |           1.98276e-05 | True          |
| preT.ETP.Th                     |   5.19913e-08 |           2.82895e-05 | True          |
| T.8SP24-.Th                     |   4.95349e-08 |           3.34413e-05 | True          |
| T.DP.69-.e17.Th.v2              |   4.53702e-08 |           3.40748e-05 | True          |
| proB.FrBC.FL                    |   4.98927e-08 |           3.65938e-05 | True          |
| MF.BM                           |   6.94982e-08 |           3.83549e-05 | True          |
| SC.MDP.BM                       |   4.9775e-08  |           3.98549e-05 | True          |
| proB.CLP.FL                     |   4.93602e-08 |           4.49107e-05 | True          |
| T.4Nve.Sp                       |   4.0749e-08  |           4.70371e-05 | True          |
| CD4.48h.LN                      |   5.13317e-08 |           4.90148e-05 | True          |
| Tgd.vg2+24ahi.e17.Th            |   3.98931e-08 |           5.05495e-05 | True          |
| T.4+8int.Th                     |   3.81905e-08 |           5.49732e-05 | True          |
| T.DP69+.Th.v2                   |   5.14471e-08 |           5.7412e-05  | True          |
| NKT.44-NK1.1-.Th                |   3.6643e-08  |           6.27774e-05 | True          |
| T.4SP69+.Th                     |   4.45522e-08 |           6.57873e-05 | True          |


These results are somewhat surprising: all categories of immune-cell are tagged as significant.  Moreover, T-helper cells are highly significant.  My guess is that this again reflects shared features of cells descended from haematopoietic stem cells.



### Corces et al. ATAC-seq data


Looking at the Corces et al. ATAC-seq dataset, we again see broad cross cell-type significance.

| Name    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------|--------------:|----------------------:|:--------------|
| Erythro |   2.20621e-06 |           5.0373e-11  | True          |
| MEP     |   5.39479e-07 |           4.68085e-08 | True          |
| CMP     |   3.95866e-07 |           7.26445e-05 | True          |
| HSC     |   3.3574e-07  |           0.000577361 | True          |
| MPP     |   2.90813e-07 |           0.00059752  | True          |
| Bcell   |   3.08234e-07 |           0.00107073  | True          |
| NK      |   3.10655e-07 |           0.00186229  | True          |
| GMP     |   2.39349e-07 |           0.00305697  | True          |


However, here we also note that erythrocyte cells are by far the most significant, followed by their direct ancestors the Megakaryocytic-Erythroid Progenitors (MEP).  This is consistent with the idea that erythrocytes are the cell type truly driving the phenotype, and the other cell types are significant to the extent that they share regulatory pathways with erythrocytes.


### Cahoy and GTEx-Brain data



| Name            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------|--------------:|----------------------:|:--------------|
| Astrocyte       |   8.37248e-09 |              0.137796 | False         |
| Oligodendrocyte |   9.32556e-10 |              0.452566 | False         |
| Neuron          |  -4.38565e-10 |              0.526099 | False         |⏎        


| Name                                    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------------------------------|--------------:|----------------------:|:--------------|
| Brain_Anterior_cingulate_cortex_(BA24)  |   1.51946e-08 |            0.00662582 | False         |
| Brain_Cerebellar_Hemisphere             |   1.18609e-08 |            0.0269974  | False         |
| Brain_Putamen_(basal_ganglia)           |   1.1835e-08  |            0.0664447  | False         |
| Brain_Cerebellum                        |   7.09048e-09 |            0.164615   | False         |
| Brain_Nucleus_accumbens_(basal_ganglia) |   5.70585e-09 |            0.170502   | False         |
| Brain_Cortex                            |   6.52218e-09 |            0.180797   | False         |
| Brain_Caudate_(basal_ganglia)           |   4.90205e-09 |            0.237016   | False         |
| Brain_Amygdala                          |   3.25693e-09 |            0.272557   | False         |



None of the neurological cell types in these datasets are significant. This is consistent with RBC volume being a non-neurological phenotype.


## Reproducing

To reproduce these results, run the analysis script [here][mecfs_bio.analysis.red_blood_cell_volume].


