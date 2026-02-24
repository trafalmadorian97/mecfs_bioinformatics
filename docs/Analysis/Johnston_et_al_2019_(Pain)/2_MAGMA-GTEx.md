---
hide:
- navigation
- toc
---

# MAGMA GTEx


Johnston et al.[@johnston2019genome] tested [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md) tissue enrichment in their study of Multisite Chronic Pain (MCP) using bulk RNAseq expression data for 53 tissues from the [GTEx database](../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex]. The significant tissues from that study, starting from most significant, were: 

1. Brain_Cortex
2. Frontal Cortex BA9 Brain_Frontal_Cortex_BA9
3. Brain_Cerebellar_Hemisphere
4. Brain_Cerebellum
5. Brain_Anterior_cingulate_cortex_BA24
6. Brain_Nucleus_accumbens_basal_ganglia
7. Brain_Caudate_basal_ganglia
8. Brain_Putamen_basal_ganglia

MAGMA gene-property analysis was repeated here using the `mecfs_bioinformatics` repo to compare the results to those in the paper. The analysis results are plotted below.

![magma-mcp-gtex](https://github.com/user-attachments/assets/273991a2-483b-4583-88f3-4314b22f97c7)

As in the original study, all of the significant tissues are brain regions. Interestingly, two additional tissues were significant in this analysis, amygdala and hippocampus, and the order of the significant tissues has changed. These differences might be attributed to the usage of different linkage disequilibrium (LD) panels as reference data for MAGMA. The original study likely used an LD panel based on the UK Biobank, as the study cohort was sampled from the UK Biobank, while the present analysis used an LD panel for Europeans from the 1000 Genomes project.
