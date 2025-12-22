# S-LDSC Analysis of DecodeME

[Stratified Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) (S-LDSC) was applied to summary statistics from GWAS 1 of [DecodeME](../../Data_Sources/MECFS/DecodeME.md).

## Reference Data Sources

I used the standard reference datasets prepared by the authors of the S-LDSC method.

- [The GTEx Project](../../Data_Sources/GTEx_Project/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../Data_Sources/Roadmap_Epigenetic_Project/Roadmap.md)
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../Data_Sources/Corces_Et_Al_Blood_Cell_Epigenetics/Corces_et_al.md).
- The [ImmGen](../../Data_Sources/Immgen_Project/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset


## Results

### GTEx and Franke lab tissue expression data
The plot and table below show the results of the application of S-LDSC to DecodeME using the GTEx and Franke lab gene expression datasets.

![s-ldsc-decodme-gene-expression](https://github.com/user-attachments/assets/77dfefc0-2717-4e45-9b1d-4f09d67b7dc5)

| Name                                                           |   Coefficient |   Coefficient_P_value | Reject Null   |
|:---------------------------------------------------------------|--------------:|----------------------:|:--------------|
| A08.186.211.Brain                                              |   4.72676e-09 |           2.47819e-07 | True          |
| A08.186.211.730.885.287.500.571.735.Visual.Cortex              |   4.53115e-09 |           1.56363e-06 | True          |
| A08.186.211.464.405.Hippocampus                                |   4.49163e-09 |           2.19491e-06 | True          |
| Brain_Cortex                                                   |   3.66873e-09 |           2.6085e-06  | True          |
| A08.186.211.730.885.287.500.270.Frontal.Lobe                   |   4.63607e-09 |           3.13473e-06 | True          |
| Brain_Nucleus_accumbens_(basal_ganglia)                        |   3.59733e-09 |           3.36967e-06 | True          |
| Brain_Frontal_Cortex_(BA9)                                     |   3.69246e-09 |           3.7678e-06  | True          |
| A08.186.211.730.885.287.500.Cerebral.Cortex                    |   4.17782e-09 |           6.18484e-06 | True          |
| A08.186.211.464.710.225.Entorhinal.Cortex                      |   4.14964e-09 |           6.71693e-06 | True          |
| A08.186.211.464.Limbic.System                                  |   4.20837e-09 |           7.01006e-06 | True          |
| Brain_Substantia_nigra                                         |   3.51471e-09 |           8.8185e-06  | True          |
| Brain_Anterior_cingulate_cortex_(BA24)                         |   3.26788e-09 |           1.0814e-05  | True          |
| Brain_Amygdala                                                 |   3.19633e-09 |           2.91008e-05 | True          |
| Brain_Hippocampus                                              |   3.25752e-09 |           3.53736e-05 | True          |
| Brain_Caudate_(basal_ganglia)                                  |   3.10888e-09 |           7.75616e-05 | True          |
| Brain_Putamen_(basal_ganglia)                                  |   3.10587e-09 |           0.000104439 | True          |
| A08.186.211.730.317.357.352.435.Hypothalamo.Hypophyseal.System |   3.61124e-09 |           0.000166442 | True          |
| A08.186.211.730.317.357.Hypothalamus                           |   3.24167e-09 |           0.00023368  | True          |
| A08.186.211.730.885.287.500.670.Parietal.Lobe                  |   3.4577e-09  |           0.000240607 | True          |
| Brain_Hypothalamus                                             |   2.84899e-09 |           0.000356557 | True          |
| A08.186.211.730.317.Diencephalon                               |   2.9906e-09  |           0.000463437 | True          |
| Brain_Cerebellar_Hemisphere                                    |   2.66524e-09 |           0.000644924 | True          |
| A08.186.211.730.885.287.249.Basal.Ganglia                      |   2.69897e-09 |           0.0010573   | True          |
| A08.186.211.865.428.Metencephalon                              |   2.44467e-09 |           0.0019363   | False         |
| A08.186.211.730.885.287.249.487.Corpus.Striatum                |   2.2671e-09  |           0.00453931  | False         |
| A08.186.211.653.Mesencephalon                                  |   2.44809e-09 |           0.00517217  | False         |
| Brain_Cerebellum                                               |   2.21062e-09 |           0.00541441  | False         |
| Brain_Spinal_cord_(cervical_c-1)                               |   2.08238e-09 |           0.00603588  | False         |


As we saw in our earlier [MAGMA analysis](MAGMA_DecodeME_Analysis.md) using the GTEx dataset, the significant tissues are all CNS-related.


### Roadmap Chromatin data

I next applied S-LDSC using the reference dataset derived from the Roadmap epigenetic project.  The results are in the plot and table below:

![s-ldsc-roadmap-chromatin](https://github.com/user-attachments/assets/79ffff5d-1693-41f5-b3b1-a54df7bd6294)


| Name                                                                     |   Coefficient |   Coefficient_P_value | Reject Null   |
|:-------------------------------------------------------------------------|--------------:|----------------------:|:--------------|
| Fetal_Brain_Female__DNase                                                |   1.02689e-07 |           4.68419e-10 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K27ac                            |   5.08217e-08 |           1.23002e-09 | True          |
| Fetal_Brain_Male__DNase                                                  |   9.55596e-08 |           2.73287e-09 | True          |
| Fetal_Brain_Male__H3K4me1                                                |   3.78077e-08 |           5.22621e-09 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K4me3                            |   1.70485e-07 |           5.98839e-09 | True          |
| Fetal_Brain_Female__H3K4me1                                              |   6.20252e-08 |           1.37543e-08 | True          |
| Brain_Anterior_Caudate__H3K27ac                                          |   3.86611e-08 |           3.06201e-08 | True          |
| Brain_Inferior_Temporal_Lobe__H3K27ac                                    |   3.30624e-08 |           2.87328e-07 | True          |
| Brain_Cingulate_Gyrus__H3K9ac                                            |   8.70527e-08 |           3.38485e-07 | True          |
| Brain_Angular_Gyrus__H3K9ac                                              |   8.39789e-08 |           9.75539e-07 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K4me1                            |   5.21768e-08 |           1.72804e-06 | True          |
| Brain_Germinal_Matrix__H3K4me3                                           |   1.26424e-07 |           1.85932e-06 | True          |
| Brain_Inferior_Temporal_Lobe__H3K4me3                                    |   1.03363e-07 |           3.9248e-06  | True          |
| Brain_Anterior_Caudate__H3K4me3                                          |   1.02781e-07 |           4.68e-06    | True          |
| Fetal_Brain_Female__H3K4me3                                              |   1.26751e-07 |           5.12867e-06 | True          |
| Brain_Angular_Gyrus__H3K4me3                                             |   1.29751e-07 |           9.79109e-06 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K9ac                             |   9.32505e-08 |           9.83818e-06 | True          |
| Brain_Cingulate_Gyrus__H3K4me3                                           |   1.06995e-07 |           1.02402e-05 | True          |
| Brain_Angular_Gyrus__H3K4me1                                             |   4.08889e-08 |           1.39444e-05 | True          |
| Brain_Angular_Gyrus__H3K27ac                                             |   3.26885e-08 |           1.89114e-05 | True          |
| Brain_Anterior_Caudate__H3K4me1                                          |   3.54259e-08 |           3.15633e-05 | True          |
| Brain_Anterior_Caudate__H3K9ac                                           |   6.35599e-08 |           3.21619e-05 | True          |
| Brain_Cingulate_Gyrus__H3K27ac                                           |   2.92786e-08 |           3.245e-05   | True          |
| Brain_Cingulate_Gyrus__H3K4me1                                           |   3.2249e-08  |           4.08821e-05 | True          |
| Cortex_derived_primary_cultured_neurospheres__H3K4me3                    |   1.30409e-07 |           5.87376e-05 | True          |
| Ganglion_Eminence_derived_primary_cultured_neurospheres__H3K4me3         |   9.46707e-08 |           8.30653e-05 | True          |
| Brain_Inferior_Temporal_Lobe__H3K9ac                                     |   5.55865e-08 |           8.69272e-05 | True          |
| Brain_Hippocampus_Middle__H3K4me3                                        |   6.50022e-08 |           0.000226226 | True          |
paiforsyth@PeterPowerDell ~/s/t/biostatistics (main)> head -n 40 assets/base_asset_store/gwas/ME_CFS/DecodeME/analysis/decode_me_gwas_1_multi_ti
ssue_chromatin_s_ldsc_cell_analysis_md_table.md
| Name                                                                     |   Coefficient |   Coefficient_P_value | Reject Null   |
|:-------------------------------------------------------------------------|--------------:|----------------------:|:--------------|
| Fetal_Brain_Female__DNase                                                |   1.02689e-07 |           4.68419e-10 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K27ac                            |   5.08217e-08 |           1.23002e-09 | True          |
| Fetal_Brain_Male__DNase                                                  |   9.55596e-08 |           2.73287e-09 | True          |
| Fetal_Brain_Male__H3K4me1                                                |   3.78077e-08 |           5.22621e-09 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K4me3                            |   1.70485e-07 |           5.98839e-09 | True          |
| Fetal_Brain_Female__H3K4me1                                              |   6.20252e-08 |           1.37543e-08 | True          |
| Brain_Anterior_Caudate__H3K27ac                                          |   3.86611e-08 |           3.06201e-08 | True          |
| Brain_Inferior_Temporal_Lobe__H3K27ac                                    |   3.30624e-08 |           2.87328e-07 | True          |
| Brain_Cingulate_Gyrus__H3K9ac                                            |   8.70527e-08 |           3.38485e-07 | True          |
| Brain_Angular_Gyrus__H3K9ac                                              |   8.39789e-08 |           9.75539e-07 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K4me1                            |   5.21768e-08 |           1.72804e-06 | True          |
| Brain_Germinal_Matrix__H3K4me3                                           |   1.26424e-07 |           1.85932e-06 | True          |
| Brain_Inferior_Temporal_Lobe__H3K4me3                                    |   1.03363e-07 |           3.9248e-06  | True          |
| Brain_Anterior_Caudate__H3K4me3                                          |   1.02781e-07 |           4.68e-06    | True          |
| Fetal_Brain_Female__H3K4me3                                              |   1.26751e-07 |           5.12867e-06 | True          |
| Brain_Angular_Gyrus__H3K4me3                                             |   1.29751e-07 |           9.79109e-06 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K9ac                             |   9.32505e-08 |           9.83818e-06 | True          |
| Brain_Cingulate_Gyrus__H3K4me3                                           |   1.06995e-07 |           1.02402e-05 | True          |
| Brain_Angular_Gyrus__H3K4me1                                             |   4.08889e-08 |           1.39444e-05 | True          |
| Brain_Angular_Gyrus__H3K27ac                                             |   3.26885e-08 |           1.89114e-05 | True          |
| Brain_Anterior_Caudate__H3K4me1                                          |   3.54259e-08 |           3.15633e-05 | True          |
| Brain_Anterior_Caudate__H3K9ac                                           |   6.35599e-08 |           3.21619e-05 | True          |
| Brain_Cingulate_Gyrus__H3K27ac                                           |   2.92786e-08 |           3.245e-05   | True          |
| Brain_Cingulate_Gyrus__H3K4me1                                           |   3.2249e-08  |           4.08821e-05 | True          |
| Cortex_derived_primary_cultured_neurospheres__H3K4me3                    |   1.30409e-07 |           5.87376e-05 | True          |
| Ganglion_Eminence_derived_primary_cultured_neurospheres__H3K4me3         |   9.46707e-08 |           8.30653e-05 | True          |
| Brain_Inferior_Temporal_Lobe__H3K9ac                                     |   5.55865e-08 |           8.69272e-05 | True          |
| Brain_Hippocampus_Middle__H3K4me3                                        |   6.50022e-08 |           0.000226226 | True          |
| Brain_Inferior_Temporal_Lobe__H3K4me1                                    |   3.38614e-08 |           0.00032817  | True          |
| Brain_Hippocampus_Middle__H3K4me1                                        |   1.98882e-08 |           0.000355901 | True          |
| Fetal_Brain_Female__H3K36me3                                             |   3.17864e-08 |           0.000550213 | True          |
| Brain_Germinal_Matrix__H3K4me1                                           |   4.83297e-08 |           0.000608534 | True          |
| skeletal_muscle_ENTEX__H3K4me1                                           |   1.17311e-08 |           0.000895252 | False         |
| Brain_Substantia_Nigra__H3K4me3                                          |   7.29864e-08 |           0.00121435  | False         |
| Placenta_Amnion__H3K36me3                                                |   2.5235e-08  |           0.00158191  | False         |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K36me3                           |   3.5113e-08  |           0.00163492  | False         |
| Brain_Hippocampus_Middle__H3K27ac                                        |   1.91901e-08 |           0.00173697  | False         |
| Ganglion_Eminence_derived_primary_cultured_neurospheres__H3K4me1         |   3.16805e-08 |           0.00215425  | False         |

Again, the strongest, and most significant associations are all with CNS cell-types

### Immgen data

Next, I applied S-LDSC using reference data from the Immgen project.

There were no significant cell types.

The cell types with the lowest p-values are shown in the table below:

| Name                            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------------------------------|--------------:|----------------------:|:--------------|
| DC.8-4-11b+.MLN                 |   2.98171e-09 |            0.00459158 | False         |
| T.4.Pa.BDC                      |   2.5939e-09  |            0.00506716 | False         |
| T.8Mem.LN                       |   2.85485e-09 |            0.00730256 | False         |
| DC.8-4-11b-.SLN                 |   2.65688e-09 |            0.0116649  | False         |
| preT.DN2-3.Th                   |   2.41341e-09 |            0.0121964  | False         |
| T.4int8+.Th                     |   2.53226e-09 |            0.012314   | False         |
| LN.TR.14w.B6                    |   2.52753e-09 |            0.0135406  | False         |
| preT.ETP-2A.Th                  |   2.50404e-09 |            0.0172496  | False         |
| MF.103-11b+.Lu                  |   2.36504e-09 |            0.0197861  | False         |
| MF.Microglia.CNS                |   2.40608e-09 |            0.0224495  | False         |
| B.Mem.Sp.v2                     |   2.06099e-09 |            0.026437   | False         |
| Tgd.vg5-.act.IEL                |   2.08167e-09 |            0.0328614  | False         |
| B.MZ.Sp                         |   2.08844e-09 |            0.0348465  | False         |
| T.8Eff.Sp.OT1.d15.LisOva        |   1.80262e-09 |            0.0357771  | False         |
| T.4SP69+.Th                     |   1.91234e-09 |            0.0394321  | False         |


### Corces et al. ATAC-seq data

The results of applying S-LDSC using the epigenetic reference data from Corces et al. ATAC-seq analysis of blood cells are shown below.  There are no significant cell types:

| Name    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------|--------------:|----------------------:|:--------------|
| Erythro |   1.83516e-08 |              0.154592 | False         |
| Bcell   |   4.60401e-09 |              0.351883 | False         |
| CLP     |   4.52542e-09 |              0.381204 | False         |
| Mono    |  -1.51999e-09 |              0.530117 | False         |
| MEP     |  -1.5488e-09  |              0.566514 | False         |
| CD4     |  -4.57403e-09 |              0.628276 | False         |
| CD8     |  -6.41332e-09 |              0.669675 | False         |
| LMPP    |  -9.21299e-09 |              0.724051 | False         |


### Cahoy and GTEx-Brain data

The next two reference datasets pertain to the nervous system.

Surprisingly, when we analyze the DecodeME results using the Cahoy dataset, the neuron cell type is not even close to being significant.  This is discordant with some of the results above, in which many CNS-related cell and tissue types were marked as significant.  Moreover, the oligodendrocyte cell type is closer to being significant than the neuron cell type.

| Name            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------|--------------:|----------------------:|:--------------|
| Oligodendrocyte |   2.12389e-09 |             0.0294642 | False         |
| Neuron          |   1.20321e-09 |             0.104572  | False         |
| Astrocyte       |  -9.45929e-10 |             0.825497  | False         |⏎ 


When we apply the S-LDSC using the GTEx brain dataset, we find the the cortex tissue type is significant:

| Name                                    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------------------------------|--------------:|----------------------:|:--------------|
| Brain_Cortex                            |   3.27938e-09 |           0.000106241 | True          |
| Brain_Frontal_Cortex_(BA9)              |   2.09994e-09 |           0.00675003  | False         |
| Brain_Anterior_cingulate_cortex_(BA24)  |   1.50154e-09 |           0.0256637   | False         |
| Brain_Nucleus_accumbens_(basal_ganglia) |   1.14427e-09 |           0.112026    | False         |
| Brain_Cerebellum                        |   8.97215e-10 |           0.182957    | False         |
| Brain_Hippocampus                       |   3.3308e-10  |           0.341041    | False         |
| Brain_Putamen_(basal_ganglia)           |   3.2917e-10  |           0.359427    | False         |
| Brain_Cerebellar_Hemisphere             |   2.45686e-10 |           0.39188     | False         


# Reproducing Analysis

To reproduce, run the [DecodeME Analysis Script][mecfs_bio.analysis.decode_me_initial_analysis]

