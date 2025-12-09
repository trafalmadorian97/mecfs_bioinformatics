# S-LDSC LDL Analysis

I applied [Stratified Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) (S-LDSC) to summary statistics from a GWAS of LDL from the [Million Veterans Program](../../Data_Sources/Million_Veterans_Program/Million_Veterans_Program.md) to identify possible key tissue and cell types affecting LDL levels.

## Reference Data Sources

I used the reference datasets recommended by the authors of the S-LSDC method.  These reference datasets are ultimately drawn from the following data sources:
- [The GTEx Project](../../Data_Sources/GTEx_Project/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../Data_Sources/Roadmap_Epigenetic_Project/Roadmap.md)
- The Cahoy Mouse Central Nervous System Dataset
- The ImmGen Project

## Results


### GTEx and Franke lab tissue expression data

When S-LSDC is applied to GWAS summary statistics using a reference dataset of cell types, S-LSDC will return cell-type $\tau_i$ coefficients together with associated $p$-values.  A large coefficient and a small $p$ value for a given cell type $i$ suggests that genes related to cell-type $i$ are over-represented in the heritability of the phenotype of interest.

The graph below shows the coefficient p-values for the cell types in the GTEx/Franke Lab dataset when L-DSC is applied to the LDL GWAS.  Cell types are grouped into categories according to the same scheme used in the original S-LDSC paper[@finucane2018heritability].

![scatter_plot](https://github.com/user-attachments/assets/c160793f-9c28-4ae8-b078-7ec1c228fe54)

Of note: the the main 3 statistically significant cell types are all liver cell types.  [Consistently with the results of MAGMA](MAGMA_LDL_Analysis.md), this suggests that the liver is  central the physiological process determining LDL levels.  Note, however, that since both S-LDSC and MAGMA use the GTEx dataset, this not truly and independent piece of evidence.

There is also one significant kidney cell type (on the right side of the graph).  It is unclear to me whether this reflects real biology, or is an artifact.

It is also of interest that several adipose-related cell types are just below the threshold of signifiance, suggesting that a larger study might implicate adipose tissue in the LDL disease process. This is consistent with the known association of high BMI with LDL levels.


Here are the top associations in tabular form:

| Name                                                           |   Coefficient |   Coefficient_P_value | Reject Null   |
|:---------------------------------------------------------------|--------------:|----------------------:|:--------------|
| A03.620.Liver                                                  |   1.5215e-08  |           0.000211744 | True          |
| A11.436.348.Hepatocytes                                        |   2.30756e-08 |           0.000254927 | True          |
| Liver                                                          |   1.82315e-08 |           0.000335455 | True          |
| Kidney_Cortex                                                  |   9.50672e-09 |           0.00045596  | True          |
| A10.165.114.830.500.750.Subcutaneous.Fat..Abdominal            |   1.13522e-08 |           0.00205337  | False         |
| A10.165.114.830.750.Subcutaneous.Fat                           |   1.11345e-08 |           0.00237814  | False         |
| A11.497.497.600.Oocytes                                        |   7.96733e-09 |           0.00312142  | False         |
| A06.407.071.Adrenal.Glands                                     |   8.71681e-09 |           0.00736218  | False         |
| A15.382.490.315.583.Neutrophils                                |   8.68251e-09 |           0.00957292  | False         |
| Pancreas                                                       |   8.30563e-09 |           0.00992878  | False         |
| A05.360.490.Germ.Cells                                         |   6.45252e-09 |           0.0110824   | False         |
| A06.407.071.140.Adrenal.Cortex                                 |   7.33265e-09 |           0.013978    | False         |
| A03.556.249.124.Ileum                                          |   1.43982e-08 |           0.0146156   | False         |


### Chromatin data