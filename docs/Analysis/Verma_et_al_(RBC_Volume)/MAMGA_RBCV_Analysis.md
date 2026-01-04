# MAGMA Analysis of RBC Volume



I applied the MAGMA pipeline to summary statistics from a GWAS of RBC from the [Million Veterans Program](../../Data_Sources/Million_Veterans_Program/Million_Veterans_Program.md).  In the gene set analysis step, I incorporated [tissue-specific RNAseq data from GTEx](../../Data_Sources/GTEx_Project/GTEx_RNAseq_Data.md) to try to link genes associated with LDL to specific tissues.


## Gene set/ Gene property analysis

The graph below shows the significant cells/tissues when MAGMA gene property analysis applied to the RBC volume GWAS using reference data from the GTEx project :

![rbc-bar-magma](https://github.com/user-attachments/assets/9d6873e4-7c44-48b7-b19d-b6e37132d4fa)

The presence of "whole blood" is intuitive.  The other four significant tissues are EBV-transformed lymphocytes, the spleen, and "Lymphode aggregate" from the small intestine.  The significance of these white-blood-cell-related tissues might be explained by the fact that erythrocytes and leukocytes share a common lineage, descending from haematopoietic stem cells (see Figure below).  Thus it could be that MAGMA with GTEx reference data is picking up genes that play a role all myeloid cells as explanatory for red-blood cell volume.


![haemopoetic-cells](https://github.com/user-attachments/assets/ccb67146-5c88-4168-be81-769a55e6b844)
Figure: By A. Rad and M. Häggström. CC-BY-SA 3.0 license

