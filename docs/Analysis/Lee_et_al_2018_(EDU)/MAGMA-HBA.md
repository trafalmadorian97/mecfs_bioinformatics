---
hide:
- navigation
- toc
---
# MAGMA HBA Analysis


I applied MAGMA to the educational attainment GWAS of Lee et al.[@lee2018gene] using scRNAseq data from the [Human Brain Atlas](../../Data_Sources/HBA_scRNAseq.md) as a reference.





In the plot below, the x-axis corresponds to HBA cluster number[@siletti2023transcriptomic], while the y-axis corresponds to the $-\log_{10}(p)$ score produced by MAGMA.  Clusters are colored according to their HBA supercluster. The dotted line denotes the Bonferroni significance threshold.  I used a conditional analysis approach based on the one described in Wanatabe et al.[@watanabe2019genetic] to identify independent clusters.  These independent clusters are labeled in plot.  


![edu-hba-magma](https://github.com/user-attachments/assets/d265ebd4-b5db-4578-a3b3-2cc02bf5f2e4)


The independent clusters are listed in more detail in the table below:

| Retained_clusters   |          P | Supercluster                   | Class auto-annotation   | Neurotransmitter auto-annotation   | Neuropeptide auto-annotation                                                     | Subtype auto-annotation   | Transferred MTG Label   | Top three regions                                                  | Top Enriched Genes                                                                                             |
|:--------------------|-----------:|:-------------------------------|:------------------------|:-----------------------------------|:---------------------------------------------------------------------------------|:--------------------------|:------------------------|:-------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------|
| Cluster234          | 1.6139e-14 | Eccentric medium spiny neuron  | NEUR                    | GABA                               | CCK CHGA CHGB NAMPT NUCB PENK SCG UBL VGF proSAAS                                | MSN-D1                    |                         | Amygdala: 75.9%, Cerebral cortex: 14.6%, Thalamus: 5.4%            | NPFFR2, ZNF736P9Y, GABRQ, LMNTD1, EYA2, AC012078.2, AC087516.2, PCDH11Y, LINC00354, NMBR                       |
| Cluster146          | 1.0347e-12 | Deep-layer intratelencephalic  | NEUR                    | VGLUT1 VGLUT2                      | ADCYAP CBLN CCK CHGA CHGB NAMPT NUCB SCG UBL VGF proSAAS                         | 0                         |                         | Amygdala: 65.4%, Cerebral cortex: 33.1%, Hypothalamus: 0.8%        | CCN4, STPG2-AS1, AC105450.1, LINC00683, ATP10A, SNTG2, UCMA, TMEM244, FGF10, LINC02008                         |
| Cluster280          | 7.3007e-09 | CGE interneuron                | NEUR                    | GABA                               | CCK CHGA CHGB CRH IGF NAMPT NUCB PNOC SCG UBL VIP proSAAS                        | INT-VIP                   |                         | Cerebral cortex: 73.7%, Amygdala: 13.3%, Hippocampus: 7.8%         | AC006305.1, LINC01416, AL391832.4, LINC01905, SP8, DCN, CHRNB3, CXCL14, CHRNA7, CRH                            |
| Cluster202          | 4.4341e-08 | Hippocampal dentate gyrus      | NEUR                    | VGLUT1                             | CBLN CHGB GRP NAMPT NUCB UBL proSAAS                                             | DG-GRAN                   |                         | Hippocampus: 100.0%, Cerebral cortex: 0.0%, Amygdala: 0.0%         | TTR, FOLR1, LINC01697, AC097528.1, AC097487.1, AL672167.1, NTF3, AL591518.1, TARID, AJ006995.1                 |
| Cluster133          | 1.1749e-07 | Upper-layer intratelencephalic | NEUR                    | VGLUT1 VGLUT2                      | ADCYAP CBLN CCK CHGA CHGB NAMPT NUCB SCG UBL VGF proSAAS                         | 0                         |                         | Cerebral cortex: 100.0%, Basal forebrain: 0.0%, Hypothalamus: 0.0% | LINC02388, AC112770.1, AC092448.1, AC126182.3, VIPR2, TRPC3, AC020637.1, IGFBP2, SLC38A11, LINC02296           |
| Cluster404          | 2.9387e-07 | Miscellaneous                  | NEUR                    | VGLUT1                             | CBLN CCK CHGA CHGB IGF NAMPT NUCB NXPH SCG SST TAC VGF proSAAS                   | 0                         |                         | Cerebral cortex: 94.4%, Amygdala: 2.5%, Basal forebrain: 1.2%      | AC016687.2, AC021134.1, LINC02378, LINC02263, LINC01915, CXCL14, AL033539.2, AC004862.1, AC073578.2, LINC02306 |
| Cluster439          | 6.1461e-07 | Midbrain-derived inhibitory    | NEUR                    | GABA                               | CBLN CHGA CHGB IGF NAMPT NPPC NUCB NXPH SCG UBL proSAAS                          | 0                         |                         | Midbrain: 52.8%, Thalamus: 31.3%, Pons: 10.7%                      | DMBX1, OTX2-AS1, OTX2, LINC01210, GATA3, AC008163.1, AL161757.4, LINC02055, SOX14, PRSS12                      |
