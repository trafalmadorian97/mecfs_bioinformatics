---
hide:
- navigation
- toc
---
# S-LDSC EDU Analysis
I applied [Stratified Linkage Disequilibrium Score Regression](../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) (S-LDSC) to summary statistics from a GWAS of educational attainment by Lee et al.[@lee2018gene]. The goal was to identify possible key tissue and cell types affecting educational attainment levels, and in so doing learn about the strengths and weaknesses of the S-LDSC method.


## Reference Data Sources

I used the standard reference datasets recommended and preprocessed by the authors of the S-LDSC method. 

- [The GTEx Project](../../Data_Sources/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../Data_Sources/Roadmap.md)
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../Data_Sources/Corces_et_al.md).
- The [ImmGen](../../Data_Sources/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset

## Results

### GTEx and Franke lab tissue expression data




The graph below shows the coefficient p values for the cell types in the GTEx/Franke Lab dataset when S-LDSC is applied to the EDU GWAS.  Cell types are grouped into categories according to the same scheme used in the original S-LDSC paper[@finucane2018heritability].

![edu_scatter](https://github.com/user-attachments/assets/a87db2e0-1d9f-4200-994f-cd1a6a3b0c7b)

We can see that at a false discovery rate of 0.01, the significant tissues are exclusively in the central nervous system.


Here is a more detailed look at the top cells/tissues:

| Name                                                           |   Coefficient |   Coefficient_P_value | Reject Null   |
|:---------------------------------------------------------------|--------------:|----------------------:|:--------------|
| A08.186.211.730.885.287.500.Cerebral.Cortex                    |   1.7625e-08  |           1.58411e-15 | True          |
| A08.186.211.464.405.Hippocampus                                |   1.62863e-08 |           4.5446e-15  | True          |
| A08.186.211.464.710.225.Entorhinal.Cortex                      |   1.68218e-08 |           6.89181e-15 | True          |
| A08.186.211.464.Limbic.System                                  |   1.67997e-08 |           3.90973e-14 | True          |
| A08.186.211.Brain                                              |   1.5237e-08  |           1.38319e-11 | True          |
| Brain_Frontal_Cortex_(BA9)                                     |   1.56443e-08 |           1.45016e-11 | True          |
| A08.186.211.730.885.287.500.571.735.Visual.Cortex              |   1.42092e-08 |           4.7507e-11  | True          |
| A08.186.211.730.885.287.500.270.Frontal.Lobe                   |   1.51812e-08 |           1.03359e-10 | True          |
| Brain_Cortex                                                   |   1.5203e-08  |           2.22635e-10 | True          |
| Brain_Anterior_cingulate_cortex_(BA24)                         |   1.3785e-08  |           4.33289e-09 | True          |
| A08.186.211.132.810.428.200.Cerebellum                         |   1.14043e-08 |           1.34142e-08 | True          |
| A08.186.211.730.885.287.500.670.Parietal.Lobe                  |   1.26008e-08 |           1.82969e-08 | True          |
| A08.186.211.865.428.Metencephalon                              |   1.10862e-08 |           5.21843e-08 | True          |
| Brain_Nucleus_accumbens_(basal_ganglia)                        |   1.27456e-08 |           5.76169e-08 | True          |
| A08.186.211.132.Brain.Stem                                     |   1.15633e-08 |           6.33204e-08 | True          |
| Brain_Hippocampus                                              |   1.20156e-08 |           2.71743e-07 | True          |
| Brain_Amygdala                                                 |   1.18184e-08 |           5.00604e-07 | True          |
| Brain_Hypothalamus                                             |   1.12805e-08 |           7.01388e-07 | True          |
| A11.872.653.Neural.Stem.Cells                                  |   1.17742e-08 |           2.32653e-06 | True          |
| Brain_Caudate_(basal_ganglia)                                  |   1.03365e-08 |           4.98279e-06 | True          |
| Brain_Putamen_(basal_ganglia)                                  |   1.04866e-08 |           6.18148e-06 | True          |
| Brain_Substantia_nigra                                         |   1.00114e-08 |           1.14614e-05 | True          |
| Brain_Cerebellum                                               |   8.55503e-09 |           1.65288e-05 | True          |
| A08.186.211.730.317.Diencephalon                               |   8.06071e-09 |           3.4515e-05  | True          |
| Brain_Cerebellar_Hemisphere                                    |   8.30636e-09 |           4.29508e-05 | True          |
| A09.371.729.Retina                                             |   8.27211e-09 |           5.85128e-05 | True          |
| A08.186.211.730.885.287.249.Basal.Ganglia                      |   8.17562e-09 |           0.000179368 | True          |
| A08.186.211.653.Mesencephalon                                  |   7.99799e-09 |           0.000264251 | True          |
| A08.186.211.730.317.357.Hypothalamus                           |   9.02518e-09 |           0.000393007 | True          |
| A08.186.211.730.317.357.352.435.Hypothalamo.Hypophyseal.System |   8.70071e-09 |           0.000577608 | True          |
| Brain_Spinal_cord_(cervical_c-1)                               |   6.92728e-09 |           0.00131103  | True          |
| A15.145.Blood                                                  |   9.14095e-09 |           0.00198242  | False         |
| A15.145.229.637.555.567.569.200.CD4.Positive.T.Lymphocytes     |   6.37535e-09 |           0.0049859   | False         |
| A11.118.637.555.567.562.440.Precursor.Cells..B.Lymphoid        |   5.41722e-09 |           0.0309822   | False         |
| A15.145.229.637.555.Leukocytes..Mononuclear                    |   4.94614e-09 |           0.0312297   | False         |
| A11.118.637.Leukocytes                                         |   4.90028e-09 |           0.0324434   | False         |
| A08.186.211.730.885.287.249.487.Corpus.Striatum                |   3.80127e-09 |           0.0324954   | False         |
| A15.378.316.580.Monocytes                                      |   5.15141e-09 |           0.035767    | False         |


Of note is that these results are not very specific: basically all CNS regions are deemed significant.  This non-specificity could be a consequence of the details of the definition of cell-type-specific genes used to create the GTEx/Franke lab datasets.  In particular, when S-LDSC is applied to the standard GTEx and Franke lab data, brain regions are not used as controls for one another.  For example, when the regression for the Cerebral Cortex is run, no other brain regions are used as controls.  This prevents us from assessing whether one brain region is better at explaining the GWAS data than another, and promotes non-specific results.


### Roadmap Chromatin data

I next applied S-LDSC to the educational attainment GWAS using reference data generated by Finucane et al.[@finucane2018heritability] from the [Roadmap Epigenetics Project](../../Data_Sources/Roadmap.md),

The following graph and table show the results:


![edu_roadmap_scatter](https://github.com/user-attachments/assets/b92bdef9-7cb0-43be-b4a0-7271893bfad6)


| Name                                                                     |   Coefficient |   Coefficient_P_value | Reject Null   |
|:-------------------------------------------------------------------------|--------------:|----------------------:|:--------------|
| Fetal_Brain_Female__H3K4me1                                              |   3.3967e-07  |           1.48993e-29 | True          |
| Fetal_Brain_Male__DNase                                                  |   4.52644e-07 |           1.46684e-25 | True          |
| Fetal_Brain_Female__DNase                                                |   4.54212e-07 |           8.02637e-25 | True          |
| Fetal_Brain_Male__H3K4me1                                                |   1.89696e-07 |           1.13027e-19 | True          |
| Fetal_Brain_Female__H3K4me3                                              |   6.47189e-07 |           2.23976e-15 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K27ac                            |   1.90671e-07 |           8.73872e-15 | True          |
| Brain_Germinal_Matrix__H3K4me3                                           |   5.93144e-07 |           1.71767e-13 | True          |
| Brain_Inferior_Temporal_Lobe__H3K27ac                                    |   1.2265e-07  |           1.02115e-12 | True          |
| Ganglion_Eminence_derived_primary_cultured_neurospheres__H3K4me3         |   4.56066e-07 |           5.04233e-12 | True          |
| Brain_Angular_Gyrus__H3K27ac                                             |   1.41502e-07 |           6.34447e-11 | True          |
| Brain_Germinal_Matrix__H3K4me1                                           |   2.78071e-07 |           6.92384e-11 | True          |
| Brain_Anterior_Caudate__H3K27ac                                          |   1.20617e-07 |           4.10239e-10 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K4me3                            |   5.32015e-07 |           4.71654e-10 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K9ac                             |   4.49235e-07 |           1.12064e-09 | True          |
| Cortex_derived_primary_cultured_neurospheres__H3K4me3                    |   5.92022e-07 |           3.65116e-09 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K4me1                            |   1.92453e-07 |           6.29012e-09 | True          |
| Fetal_Brain_Female__H3K36me3                                             |   1.35901e-07 |           1.00633e-08 | True          |
| Brain_Angular_Gyrus__H3K4me3                                             |   4.09788e-07 |           1.86578e-08 | True          |
| Brain_Anterior_Caudate__H3K4me3                                          |   3.15918e-07 |           3.9166e-08  | True          |
| Brain_Angular_Gyrus__H3K9ac                                              |   2.47642e-07 |           4.53828e-08 | True          |
| Ganglion_Eminence_derived_primary_cultured_neurospheres__H3K4me1         |   1.55274e-07 |           5.48653e-08 | True          |
| Brain_Inferior_Temporal_Lobe__H3K4me3                                    |   3.60108e-07 |           6.96504e-08 | True          |
| Brain_Anterior_Caudate__H3K4me1                                          |   1.21433e-07 |           2.15232e-07 | True          |
| Cortex_derived_primary_cultured_neurospheres__H3K4me1                    |   1.0066e-07  |           2.91758e-07 | True          |
| Brain_Anterior_Caudate__H3K9ac                                           |   2.20957e-07 |           6.50058e-07 | True          |
| Brain_Cingulate_Gyrus__H3K4me3                                           |   3.26012e-07 |           1.20709e-06 | True          |
| Brain_Germinal_Matrix__H3K36me3                                          |   1.46017e-07 |           1.23627e-06 | True          |
| Brain_Cingulate_Gyrus__H3K9ac                                            |   2.32592e-07 |           2.80626e-06 | True          |
| Brain_Inferior_Temporal_Lobe__H3K9ac                                     |   1.75209e-07 |           4.73495e-06 | True          |
| Brain_Cingulate_Gyrus__H3K27ac                                           |   8.32185e-08 |           6.05097e-06 | True          |
| Brain_Angular_Gyrus__H3K4me1                                             |   1.02775e-07 |           7.61444e-06 | True          |
| Brain_Dorsolateral_Prefrontal_Cortex__H3K36me3                           |   1.32893e-07 |           3.58725e-05 | True          |
| Brain_Cingulate_Gyrus__H3K4me1                                           |   8.40684e-08 |           4.28436e-05 | True          |
| skeletal_muscle_ENTEX__H3K4me1                                           |   3.55999e-08 |           6.17961e-05 | True          |
| Brain_Hippocampus_Middle__H3K27ac                                        |   6.5831e-08  |           6.30559e-05 | True          |
| Brain_Hippocampus_Middle__H3K4me3                                        |   1.87842e-07 |           0.000196819 | True          |
| Brain_Hippocampus_Middle__H3K4me1                                        |   5.8769e-08  |           0.000200932 | True          |
| Brain_Inferior_Temporal_Lobe__H3K4me1                                    |   9.74402e-08 |           0.000227256 | True          |
| Brain_Inferior_Temporal_Lobe__H3K36me3                                   |   9.70243e-08 |           0.000449288 | True          |
| Pancreatic_Islets__H3K27ac                                               |   1.52361e-07 |           0.00230255  | False         |
| Brain_Angular_Gyrus__H3K36me3                                            |   1.0321e-07  |           0.00275453  | False         |
| Brain_Anterior_Caudate__H3K36me3                                         |   7.73984e-08 |           0.00310657  | False         |
| Esoph-GJ_ENTEX__H3K27ac                                                  |   4.55769e-08 |           0.00313903  | False         |
| Colon-Sigm_ENTEX__H3K36me3                                               |   6.8338e-08  |           0.00367321  | False         |
| Esoph-Muscularis_ENTEX__H3K4me3                                          |   1.53155e-07 |           0.00454078  | False         |
| Foreskin_Fibroblast_Primary_Cells_skin02__H3K4me3                        |   6.39511e-08 |           0.00565107  | False         |
| Esoph-Muscularis_ENTEX__H3K27ac                                          |   4.31421e-08 |           0.00586156  | False         |
| Spleen__H3K36me3                                                         |   4.78941e-08 |           0.00599588  | False         |
| Esoph-GJ_ENTEX__H3K4me3                                                  |   1.3337e-07  |           0.00875963  | False         |
| Brain_Substantia_Nigra__H3K4me1                                          |   4.42562e-08 |           0.0114508   | False         |
| Brain_Cingulate_Gyrus__H3K36me3                                          |   6.88419e-08 |           0.0126657   | False         |
| Esoph-Mucosa_ENTEX__H3K36me3                                             |   5.47652e-08 |           0.0216056   | False         |
| Fetal_Brain_Male__H3K4me3                                                |   4.21314e-07 |           0.0237894   | False         |
| Colon-Sigm_ENTEX__H3K27ac                                                |   3.46076e-08 |           0.0241773   | False         |
| Brain_Substantia_Nigra__H3K4me3                                          |   1.31624e-07 |           0.0253406   | False         |
| Heart-LV_ENTEX__H3K4me1                                                  |   4.07981e-08 |           0.0254759   | False         |
| Nerve-Tibial_ENTEX__H3K4me3                                              |   1.51766e-07 |           0.026091    | False         |
| Brain_Hippocampus_Middle__H3K36me3                                       |   5.12254e-08 |           0.0268535   | False         |

Consistent with the GTEX/Franke results, here many regions of the brain are deemed significant.  Because the Roadmap epigenetics dataset is entirely separate from the GTEx gene expression dataset, these results can be seen as an independent piece of evidence for the significance of the brain to educational attainment.

Similarly to the GTEx/Franke case, the non-specificity of the brain results may also be a consequence of the lack of brain controls in the regressions used to define cell-type-specific regions.

The presence of a single significant skeletal muscle cell type is unexpected.  It is unclear if this is a real effect or an artifact.


### ImmGen data
The next step is to use the S-LDSC reference data derived from the ImmGen project.  There are no significant hits.  The top non-significant results are:

| Name                            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------------------------------|--------------:|----------------------:|:--------------|
| Tgd.vg2+24ahi.Th                |   8.06392e-09 |            0.00297249 | False         |
| NKT.4+.Sp                       |   7.13596e-09 |            0.00592045 | False         |
| T.8Eff.Sp.OT1.d15.LisOva        |   5.48845e-09 |            0.0103446  | False         |
| T.DPsm.Th                       |   6.48315e-09 |            0.0113907  | False         |
| B.FrE.BM                        |   5.14779e-09 |            0.0116146  | False         |
| B.T2.Sp                         |   5.60899e-09 |            0.0119636  | False         |
| CD8.1h.LN                       |   6.08251e-09 |            0.0120099  | False         |
| NK.Sp                           |   5.82532e-09 |            0.0122281  | False         |


### Corces et al. ATAC-seq data

The next step is to use the Corces et al. dataset from ATAC-seq applied to blood cells.

Again, there were no significant results.  The top non-significant results are:

| Name    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:--------|--------------:|----------------------:|:--------------|
| Erythro |   4.19012e-08 |              0.148498 | False         |
| MEP     |   2.21546e-08 |              0.285599 | False         |
| Bcell   |   1.09094e-08 |              0.360172 | False         |
| MPP     |   7.66168e-09 |              0.402965 | False         |
| CMP     |   5.85289e-09 |              0.426029 | False         |
| CLP     |   6.38071e-09 |              0.431636 | False         |
| HSC     |  -1.82334e-09 |              0.521691 | False         |
| GMP     |  -1.83842e-08 |              0.735746 | False         |


### Cahoy and GTEx-Brain data

The next two reference datasets pertain to the nervous system.


The Cahoy CNS datasets points to neurons as the key site of activity driving educational attainment, as opposed to oligodendrocytes or astrocytes:

| Name            |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------|--------------:|----------------------:|:--------------|
| Neuron          |   1.23199e-08 |           4.25819e-08 | True          |
| Oligodendrocyte |  -2.70708e-09 |           0.862452    | False         |
| Astrocyte       |  -2.84009e-09 |           0.902723    | False         |⏎


Using the GTEx-brain dataset points to the cortex as the key region of the brain for educational attainment:

| Name                                    |   Coefficient |   Coefficient_P_value | Reject Null   |
|:----------------------------------------|--------------:|----------------------:|:--------------|
| Brain_Frontal_Cortex_(BA9)              |   9.64677e-09 |           6.32016e-07 | True          |
| Brain_Cortex                            |   9.55629e-09 |           7.00572e-06 | True          |
| Brain_Anterior_cingulate_cortex_(BA24)  |   6.67343e-09 |           0.000201072 | True          |
| Brain_Cerebellar_Hemisphere             |   4.16361e-09 |           0.02365     | False         |
| Brain_Nucleus_accumbens_(basal_ganglia) |   4.34605e-09 |           0.0319243   | False         |
| Brain_Cerebellum                        |   2.81524e-09 |           0.107455    | False         |
| Brain_Putamen_(basal_ganglia)           |  -2.16168e-10 |           0.535069    | False         |
| Brain_Amygdala                          |  -1.67611e-09 |           0.828353    | False         |



The strong contrast between these results (which are highly specific to the cortex) and the results from the general GTEx dataset  [above](#gtex-and-franke-lab-tissue-expression-data) (which mark all brain regions) is explained in the Methods section of the Finucane et al. paper[@finucane2018heritability].  As explained there, in the GTEx brain S-LDSC reference dataset, a gene is associated with a particular brain region only if that gene is a strongly significant predictor in a regression  aimed at distinguishing the given brain region from all other brain regions. Thus we should expect the GTEx brain dataset to be more informative than the general GTEx dataset when we aim to determine the specific brain region driving the educational attainment phenotype.



## Reproducing Analysis

To reproduce, run the script [here][mecfs_bio.analysis.lee_educational_attainment_analysis].


