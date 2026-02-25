---
hide:
- toc
---
# MAGMA HBA Analysis
I applied MAGMA to the Diastolic Blood Pressure GWAS of Keaton et al.[@keaton2024genome]  using scRNAseq data from the [Human Brain Atlas](../../Data_Sources/HBA_scRNAseq.md) as a reference.

## Results
The results are plotted below:

![hba_magma_indep_plot](https://github.com/user-attachments/assets/751f2f5e-6d85-4c1b-8604-57d72a1dd931)

The x-axis corresponds to HBA cluster number[@siletti2023transcriptomic], while the y-axis corresponds to the $-\log_{10}(p)$ score produced by MAGMA.  Clusters are colored according to their HBA supercluster. The dotted line denotes the Bonferroni significance threshold.  I used a conditional analysis approach based on the one described in Wanatabe et al.[@watanabe2019genetic] to identify independent clusters.  These 3 independent clusters are labeled in plot.  I have also listed them in the table below, together with some cluster-annotations from Duncan et al.[@duncan2025mapping].


| Retained_clusters   |          P | Supercluster   | Class auto-annotation   |   Neurotransmitter auto-annotation | Neuropeptide auto-annotation | Subtype auto-annotation   | Transferred MTG Label   | Top three regions                                                 | Top Enriched Genes                                                                       |
|:--------------------|-----------:|:---------------|:------------------------|-----------------------------------:|:-----------------------------|:--------------------------|:------------------------|:------------------------------------------------------------------|:-----------------------------------------------------------------------------------------|
| Cluster20           | 7.291e-13  | Vascular       | VSMC                    |                                  0 | AGT NAMPT UBL                | 0                         |                         | Thalamus: 21.8%, Cerebral cortex: 13.7%, Medulla: 12.4%           |[ MYL9](https://web.archive.org/web/20250814130420/https://www.genecards.org/cgi-bin/carddisp.pl?gene=MYL9), CARMN, ADIRF, FHL5, ACTA2, VIM, TBX18, MUSTN1, TBX2, NOTCH3                        |
| Cluster57           | 4.5477e-08 | Astrocyte      | ASTRO                   |                                  0 | [AGT](https://www.uniprot.org/uniprotkb/P01019/entry)                      | 0                         | Astro                   | Hippocampus: 66.8%, Basal forebrain: 25.3%, Cerebral cortex: 4.0% | [CD38](https://en.wikipedia.org/wiki/CD38), OBI1-AS1, FGFR3, AL627316.1, AC012405.1, AL137139.2, TNC, LINC02649, GLI3, SLC14A1 |
| Cluster16           | 1.4648e-07 | Vascular       | ENDO                    |                                  0 | NAMPT NUCB                   | CAP                       | Endo                    | Cerebral cortex: 25.5%, Hippocampus: 17.8%, Midbrain: 8.4%        | [CLDN5](https://en.wikipedia.org/wiki/CLDN5), MECOM, VWF, FLT1, ABCG2, FLI1, NOSTRIN, LINC02147, TGM2, CFH                      |⏎



The significance of Cluster20 and Cluster16, which are auto-annotated as vascular smooth muscle cells (VSMC) and endothelial cells (ENDO) is unsurprising: the importance of vascular cell types to blood pressure is well-known.  This finding probably does not reflect brain-specific physiology, but rather the general causal connection between vascular cells and blood pressure.

In contrast, the significance of the astrocyte cell type does reflect brain-specific physiology.  As seen in the neuropeptide auto-annotation column in the table above, astrocytes are producers of angiotensinogen ([AGT](https://www.uniprot.org/uniprotkb/P01019/entry)), precursor to [angiotensin](https://en.wikipedia.org/wiki/Angiotensin), which in turn is a key blood-pressure-modulating hormone.  There are numerous confirmatory experimental studies demonstrating the role of astrocytes in the control of blood pressure[@stern2016astrocytes;@verma2025emerging].