---
hide:
  - navigation
  - toc
---
# MAGMA HBA Analysis
I applied MAGMA to the Liu et al.[@liu2023genetic] GWAS of inflammatory bowel disease  using scRNAseq data from the [Human Brain Atlas](../../Data_Sources/HBA_scRNAseq.md) as a reference.

## Results
The results are plotted below:



![hba-ibd-magma-plot](https://github.com/user-attachments/assets/3eab64b8-14af-4f80-92db-7c8a3631d5d6)


The x-axis corresponds to HBA cluster number[@siletti2023transcriptomic], while the y-axis corresponds to the $-\log_{10}(p)$ score produced by MAGMA.  Clusters are colored according to their HBA supercluster. The dotted line denotes the Bonferroni significance threshold.  I used a conditional analysis approach based on the one described in Wanatabe et al.[@watanabe2019genetic] to identify independent clusters.  These 2 independent clusters are labeled in plot.  I have also listed them in the table below, together with some cluster-annotations from Duncan et al.[@duncan2025mapping].


| Retained_clusters   |          P | Supercluster   | Class auto-annotation   |   Neurotransmitter auto-annotation | Neuropeptide auto-annotation   |   Subtype auto-annotation | Transferred MTG Label   | Top three regions                                    | Top Enriched Genes                                                    |
|:--------------------|-----------:|:---------------|:------------------------|-----------------------------------:|:-------------------------------|--------------------------:|:------------------------|:-----------------------------------------------------|:----------------------------------------------------------------------|
| Cluster4            | 4.6356e-14 | Microglia      | MGL                     |                                  0 | NAMPT                          |                         0 | Micro-PVM               | Spinal cord: 31.8%, Pons: 26.0%, Medulla: 13.2%      | SRGN, RGS1, GPR183, CD69, HLA-DRA, OLR1, TNFRSF1B, IFI30, CXCR4, CD74 |
| Cluster1            | 3.3757e-10 | Miscellaneous  | TCELL                   |                                  0 | 0                              |                         0 |                         | Midbrain: 15.0%, Basal forebrain: 14.0%, Pons: 13.2% | CD2, IL7R, PTPRC, SLFN12L, IL32, CCL5, GRAP2, RUNX3, CD69, CD3E       |⏎


- Cluster2 (T-cells) is consistent with the known biology of inflammatory bowel disease.  The importance of T-cells to inflammatory bowel disease is well known.
- Cluster4 (Microglia) is more surprising.  While there are some ("Gut-Brain-Axis") theories the give a role to the CNS in the development of Crohn's disease, it seems to me more likely that the significance of this cluster in MAGMA analysis is an artifact. Probably, IBD is associated with certain myeloid transcriptional programs, and microglia are the best representatives of these transcriptional programs in the human brain atlas dataset. 