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

![scatter_plot](https://github.com/user-attachments/assets/e074b62c-ac45-4ac7-b436-738679e0cc6d)

Of note: the only 3 statistically significant cell types are all liver cell types.  [Consistently with the results of MAGMA](MAGMA_LDL_Analysis.md), this suggests that the liver is  central the physiological process determining LDL levels.  Note, however, that since both S-LDSC and MAGMA use the GTEx dataset, this not truly and independent piece of evidence.

It is also of interest that several adipose-related cell types are just below the threshold of signifiance, suggesting that a larger study might implicate adipose tissue in the LDL disease process.


### Chromatin data