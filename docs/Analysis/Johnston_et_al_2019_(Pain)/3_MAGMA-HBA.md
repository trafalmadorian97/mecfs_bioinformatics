---
hide:
- toc
---
# MAGMA HBA
MAGMA was used to test the Johnston 2019[@johnston2019genome] GWAS of Multisite Chronic Pain (MCP) using scRNAseq data from the [Human Brain Atlas](../../Data_Sources/HBA_scRNAseq.md)[@siletti2023transcriptomic] as a reference.

## Results
The results are plotted below:

![hba-magma-mcp](https://github.com/user-attachments/assets/0c098a3b-7543-4c1c-9f1f-1b3440b56be7)

The x-axis corresponds to HBA cluster number, while the y-axis corresponds to the $-\log_{10}(p)$ score produced by MAGMA.  Clusters are colored according to their HBA supercluster. The dotted line denotes the Bonferroni significance threshold.

The most significant cell types are in the [Caudate ganglionic eminences (CGE) interneuron](https://www.proteinatlas.org/humanproteome/single+cell/single+nuclei+brain/neuronal+cells#cge_interneuron) supercluster.

After applying a stepwise selection algorithm similar to the one describe in Watanabe et al.[@watanabe2019genetic], three clusters were identified as independently significant:

| Retained_clusters   |          P | Supercluster                  | Class auto-annotation   | Neurotransmitter auto-annotation   | Neuropeptide auto-annotation                                                         | Subtype auto-annotation   | Transferred MTG Label   | Top three regions                                                | Top Enriched Genes                                                                                      |
|:--------------------|-----------:|:------------------------------|:------------------------|:-----------------------------------|:-------------------------------------------------------------------------------------|:--------------------------|:------------------------|:-----------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------|
| Cluster286          | 1.8883e-07 | CGE interneuron               | NEUR                    | GABA                               | CBLN CCK CHGA CHGB CORT CRH IGF NAMPT NPW NPY NUCB NXPH PNOC SCG SST UBL VGF proSAAS | INT-LAMP5                 | Lamp5                   | Cerebral cortex: 70.6%, Hippocampus: 13.8%, Amygdala: 10.2%      | NDNF, CXCL14, AC005064.1, RELN, NMBR, LINC00298, MYH11, AC092819.3, NR2E1, SLC25A48                     |
| Cluster132          | 1.1917e-05 | Miscellaneous                 | NEUR                    | VGLUT1 VGLUT2                      | CBLN CCK CHGA CHGB NAMPT NUCB SCG UBL VGF proSAAS                                    | 0                         |                         | Cerebral cortex: 98.4%, Hippocampus: 0.8%, Midbrain: 0.8%        | CD53, DPP4, AC060765.1, LINC01721, AC112770.1, AL157769.1, VIPR2, PCSK1, AC016687.2, QRFPR              |
| Cluster137          | 1.7584e-05 | Deep-layer intratelencephalic | NEUR                    | VGLUT1 VGLUT2                      | CBLN CCK CHGA CHGB CRH NAMPT NUCB SCG UBL UCN VGF proSAAS                            | 0                         | L5 IT                   | Cerebral cortex: 94.4%, Basal forebrain: 2.1%, Hippocampus: 1.9% | LINC02196, AC099517.1, AC079380.1, AC073578.2, ADGRL4, PKD2L1, AC021134.1, AL450352.1, TRIM22, ARHGAP15 |






