# MAGMA HBA Analysis

I applied MAGMA to the PGC 2022 GWAS of schizophrenia using single-cell RNAseq read-counts from the [Human Brain Atlas](../../Data_Sources/Human_Brain_Atlas/HBA_scRNAseq.md)[@siletti2023transcriptomic] as reference data. This was an attempt to verify that I could reproduce the results of Duncan et al.[@duncan2025mapping].

## Results

The results are plotted below:

![sch-human-brain-magma](https://github.com/user-attachments/assets/55268088-a7c2-446e-aec3-0aca6ee5b9ed)

These results closely match those of Duncan et al.  Like Duncan et al., we find that the most significant cell type for Schizophrenia is the Somatostatin MGE Interneuron.


When we apply a technique based on Watanabe et al.'s to identify independently significant clusters, we get the following results:


| Retained_clusters   |          P | Supercluster                  | Class auto-annotation   | Neurotransmitter auto-annotation   | Neuropeptide auto-annotation                                                     | Subtype auto-annotation   | Transferred MTG Label   | Top three regions                                                | Top Enriched Genes                                                                                             |
|:--------------------|-----------:|:------------------------------|:------------------------|:-----------------------------------|:---------------------------------------------------------------------------------|:--------------------------|:------------------------|:-----------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------|
| Cluster239          | 5.4105e-16 | MGE interneuron               | NEUR                    | GABA                               | CCK CHGA CHGB CORT IGF NAMPT NMU NPPC NUCB NXPH PNOC SCG SST TAC UBL VGF proSAAS | INT-SST                   |                         | Cerebral cortex: 95.4%, Hippocampus: 1.5%, Basal forebrain: 1.2% | AC022905.1, NMU, MAFB, SCARA5, PRDM1, TRBC2, MYO5B, PDCL2, BLM, PNOC                                           |
| Cluster280          | 1.8714e-15 | CGE interneuron               | NEUR                    | GABA                               | CCK CHGA CHGB CRH IGF NAMPT NUCB PNOC SCG UBL VIP proSAAS                        | INT-VIP                   |                         | Cerebral cortex: 73.7%, Amygdala: 13.3%, Hippocampus: 7.8%       | AC006305.1, LINC01416, AL391832.4, LINC01905, SP8, DCN, CHRNB3, CXCL14, CHRNA7, CRH                            |
| Cluster132          | 2.4393e-13 | Miscellaneous                 | NEUR                    | VGLUT1 VGLUT2                      | CBLN CCK CHGA CHGB NAMPT NUCB SCG UBL VGF proSAAS                                | 0                         |                         | Cerebral cortex: 98.4%, Hippocampus: 0.8%, Midbrain: 0.8%        | CD53, DPP4, AC060765.1, LINC01721, AC112770.1, AL157769.1, VIPR2, PCSK1, AC016687.2, QRFPR                     |
| Cluster231          | 1.9584e-12 | Eccentric medium spiny neuron | NEUR                    | GABA                               | CHGB NAMPT NUCB PDYN PENK SCG VGF proSAAS                                        | MSN-D1                    |                         | Amygdala: 60.2%, Basal forebrain: 16.4%, Cerebral cortex: 13.7%  | LMNTD1, NPFFR2, AC012078.2, ZNF736P9Y, GABRQ, EYA2, PCDH11Y, CAPSL, DRAXIN, ITGA6                              |
| Cluster404          | 6.6585e-08 | Miscellaneous                 | NEUR                    | VGLUT1                             | CBLN CCK CHGA CHGB IGF NAMPT NUCB NXPH SCG SST TAC VGF proSAAS                   | 0                         |                         | Cerebral cortex: 94.4%, Amygdala: 2.5%, Basal forebrain: 1.2%    | AC016687.2, AC021134.1, LINC02378, LINC02263, LINC01915, CXCL14, AL033539.2, AC004862.1, AC073578.2, LINC02306 |
| Cluster407          | 7.0269e-06 | Amygdala excitatory           | NEUR                    | VGLUT2 VGLUT3                      | CART CCK CHGA CHGB NAMPT NUCB NXPH SCG TAC UBL VGF proSAAS                       | 0                         |                         | Amygdala: 72.6%, Cerebral cortex: 23.5%, Hippocampus: 2.7%       | TFAP2C, PAPPA2, CARTPT, BARHL2, LINC01798, MEIS1-AS3, AC006065.4, LHX9, GALR1, POU2F3                          |
| Cluster456          | 1.6412e-05 | Thalamic excitatory           | NEUR                    | VGLUT2                             | ADCYAP CART CBLN CCK CHGA CHGB CRH NAMPT NUCB NXPH SCG UBL VGF proSAAS           | 0                         |                         | Thalamus: 99.8%, Hippocampus: 0.2%, Cerebral cortex: 0.0%        | ST8SIA6, AC012535.1, AC093334.1, ST8SIA6-AS1, AP000331.1, AC006487.1, SHOX2, ANXA3, AC087482.1, VANGL1         |



Here is a the plot of significant clusters again, this time with independently significant clusters labeled:

![sch-hba-cluster-plot](https://github.com/user-attachments/assets/f859d2e5-b507-4c2e-b658-0f89f6fdef97)

It is interesting to observe the fairly significant overlap with the significant clusters found in the [analysis of the educational attainment GWAS](../Lee_et_al_2018_(EDU)/MAGMA-HBA.md).