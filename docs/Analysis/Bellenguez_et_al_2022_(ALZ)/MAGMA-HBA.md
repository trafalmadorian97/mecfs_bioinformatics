---
hide:
- navigation
- toc
---
# MAGMA HBA
I applied MAGMA to the Alzheimer's GWAS of Bellenguez et al.[@bellenguez2022new] using scRNAseq data from the [Human Brain Atlas](../../Data_Sources/HBA_scRNAseq.md)[@siletti2023transcriptomic] as a reference.

## Results
The results are plotted below:

![hba-magma-alz](https://github.com/user-attachments/assets/a8940dff-ce73-4480-9e78-95e49d4d357d)

The x-axis corresponds to HBA cluster number, while the y-axis corresponds to the $-\log_{10}(p)$ score produced by MAGMA.  Clusters are colored according to their HBA supercluster. The dotted line denotes the Bonferroni significance threshold.


Almost all of the Bonferroni-significant cell types are in the microglia supercluster.  The only exception is a single monocyte cell type.  This is consistent with the theory of Alzheimer's as a primarily microglia-driven disease  (see pg. 1205 of Kandel et al.[@kandel2021principles]).


Applying a stepwise selection algorithm similar to the one describe in Watanabe et al.[@watanabe2019genetic] identifies only a single independently significant cluster: Cluster 12, which consists of microglial cells.