# S-LDSC


I used [Stratified Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) (S-LDSC) to analyze summary statistics from a GWAS of heart rate recovery after exercise 

I used the standard reference data sources to associate chromatin regions with cell or tissue types. These reference data sources were:

- [The GTEx Project](../../Data_Sources/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../Data_Sources/Roadmap.md)
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../Data_Sources/Corces_et_al.md).
- The [ImmGen](../../Data_Sources/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset

Surprisingly, there were no significant tissue or cell types in any of these datasets. This may indicate a poor match between the cell and tissue types present in these datasets and the physiology of heart rate recovery.