---
hide:
- toc
---
# HBA MAGMA


I applied MAGMA to the heart rate recovery GWAS of Verweiji et al.[@verweij2018genetic] using scRNAseq data from the [Human Brain Atlas](../../Data_Sources/HBA_scRNAseq.md) as a reference.





In the plot below, the x-axis corresponds to HBA cluster number[@siletti2023transcriptomic], while the y-axis corresponds to the $-\log_{10}(p)$ score produced by MAGMA.  Clusters are colored according to their HBA supercluster. The dotted line denotes the Bonferroni significance threshold.  I used a conditional analysis approach based on the one described in Wanatabe et al.[@watanabe2019genetic] to identify independent clusters.   A single independent cluster was found, which is labeled in the plot, and listed in detail the table below.




[//]: # (<iframe src="../../../_figs/verweiji_et_al_hrr_hba_magma_plot_extracted.html" style="width:100%; height:600px; border:none;"></iframe>)



![verweiji_independent_plot](https://github.com/user-attachments/assets/63dbf730-2aca-4f6d-9d19-f32de74a788d)

| Retained_clusters   |          P | Supercluster   | Class auto-annotation   | Neurotransmitter auto-annotation   | Neuropeptide auto-annotation                                               |   Subtype auto-annotation | Transferred MTG Label   | Top three regions                           | Top Enriched Genes                                                              |
|:--------------------|-----------:|:---------------|:------------------------|:-----------------------------------|:---------------------------------------------------------------------------|--------------------------:|:------------------------|:--------------------------------------------|:--------------------------------------------------------------------------------|
| Cluster398          | 1.5129e-07 | Splatter       | NEUR                    | VGLUT2                             | ADCYAP CART CBLN CCK CHGA CHGB GAL NAMPT NUCB NXPH SCG TAC UBL VGF proSAAS |                         0 |                         | Medulla: 86.6%, Pons: 12.8%, Midbrain: 0.2% | PHOX2B-AS1, PHOX2B, LINC00682, HOXB3, SLC5A7, SLC6A2, GAL, TH, PHOX2A, HOTAIRM1 |⏎   

This result is striking.  Notice that the cell type in the table above is found mostly in the [medulla oblongata](https://en.wikipedia.org/wiki/Medulla_oblongata), a region of the brainstem. The brainstem is a center of autonomic nervous system control. Given that the autonomic nervous system is key to the control of heart rate, it makes physiologic sense that genes associated with heart rate recovery after exercise should be expressed in an autonomic control center.