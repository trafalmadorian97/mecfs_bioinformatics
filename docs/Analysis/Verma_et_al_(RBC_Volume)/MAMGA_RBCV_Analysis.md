---
hide:
- toc
---
# MAGMA Analysis of RBC Volume



I applied the MAGMA pipeline to summary statistics from a GWAS of red blood cell (RBC) volume from the [Million Veterans Program](../../Data_Sources/Million_Veterans_Program.md)[@verma2024diversity].  I incorporated [tissue-specific RNAseq data from GTEx](../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex] to try to link genes associated with RBC volume to specific tissues.


## Gene set/ Gene property analysis

The graph below shows the significant cells/tissues when MAGMA gene property analysis applied to the RBC volume GWAS using reference data from the GTEx project:

![rbc-bar-magma](https://github.com/user-attachments/assets/9d6873e4-7c44-48b7-b19d-b6e37132d4fa)

The presence of "whole blood" is intuitive.  The other four significant tissues are EBV-transformed lymphocytes, the spleen, and "Lymphode aggregate" from the small intestine.  The significance of these white-blood-cell-related tissues might be explained by their similarity to erythrocytes. Erythrocytes and leukocytes share a common lineage, descending from haematopoietic stem cells (see Figure below).  Thus it could be that MAGMA is picking up genes that play a role in all myeloid cells.


![haemopoetic-cells](https://github.com/user-attachments/assets/ccb67146-5c88-4168-be81-769a55e6b844)
Figure: By A. Rad and M. Häggström. CC-BY-SA 3.0 license

