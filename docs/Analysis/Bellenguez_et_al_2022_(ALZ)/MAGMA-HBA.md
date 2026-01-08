# MAGMA HBA
I applied MAGMA to the Alzheimer's GWAS of Bellenguez et al. using scRNAseq data from the [Human Brain Atlas](../../Data_Sources/Human_Brain_Atlas/HBA_scRNAseq.md) as a reference.

## Results
The results are plotted below:

![hba-magma-alz](https://github.com/user-attachments/assets/a8940dff-ce73-4480-9e78-95e49d4d357d)

The x-axis corresponds to HBA cluster number, while the y axis corresponds to the $-\log_{10}(p)$ value score produced by magma.  Clusters are colored according to their HBA supercluster. The dotted line denotes the Bonferroni significance threshold.


Almost all of the Bonferroni-significant cell types are in the microglia supercluster.  The only except is a single monocyte cell type.  This is highly consistent with the theory of Alzheimer's as a primarily microglia-driven disease  (see pg. 1205 of Kandel et al.,[@kandel2021principles]).