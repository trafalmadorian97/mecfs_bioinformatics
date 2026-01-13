# HBA scRNAseq Dataset`

The Human Brain Atlas (HBA)[@siletti2023transcriptomic] is a recent project aiming to survey the diversity of cell types across the brain through single-cell RNA sequencing (scRNAseq).

The HBA dataset was derived from scRNAseq performed on 106 dissections of the brains of 3 post-mortem donors.  The result was gene expression profiles from over 3 million brain cells.





## Clustering

HBA analysts applied hierarchical clustering algorithms to the RNA expression levels of these ~3 million cells to create clusters at 3 scales:  subclusters, clusters, and superclusters. 

The result is that "brain cell types" in the HBA dataset are defined in a mostly data-driven way.

The Figure below is from Silleti et al.[@siletti2023transcriptomic].  The left panel illustrates the dissected regions of the brain.  The right panel shows a t-SNE plot of brain cells, with each cell colored by its supercluster.
![hba-plot](https://github.com/user-attachments/assets/146c89c8-da31-41c2-a56c-903ef3de0f28)


## Utility 

The HBA dataset can be combined with [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md) to associate GWAS phenotypes with HBA-defined brain cell types.