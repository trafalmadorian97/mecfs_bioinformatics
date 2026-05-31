---
hide:
- toc
---
# MAGMA HBA
MAGMA was used to test the UK BioBank GWAS of migraine[@uk2025whole][@GCST90473326_migraine] using scRNAseq data from the [Human Brain Atlas](../../../Data_Sources/HBA_scRNAseq.md)[@siletti2023transcriptomic] as a reference.

## Results
The results are plotted below:

{{ plotly_embed("docs/_figs/uk_biobank_2025_migraine_eur_hba_magma_plot_extracted.html", id="uk_biobank_2025_migraine_eur_hba_magma_plot_extracted", caption="The x-axis corresponds to HBA cluster number, while the y-axis corresponds to the -log10(p) score produced by MAGMA.  Clusters are colored according to their HBA supercluster. The dotted line denotes the Bonferroni significance threshold." ) }}

In line with the MAGMA tissue results, no brain cell types were significant after multiple test correction. The most significant cell type was [LAMP5-LHX6 and Chandelier cell cluster](https://www.proteinatlas.org/humanproteome/single+cell/single+nuclei+brain/neuronal+cells#lamp5_lhx6_and_chandelier).