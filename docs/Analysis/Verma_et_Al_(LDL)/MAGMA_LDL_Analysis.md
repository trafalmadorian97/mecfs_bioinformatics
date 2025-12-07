

# MAGMA Analysis of LDL



## MAGMA Results

I applied the MAGMA pipeline to summary statistics from a GWAS of LDL from the [Million Veterans Program](../../Data_Sources/Million_Veterans_Program/Million_Veterans_Program.md).  In the gene set analysis step, I incorporated [tissue-specific RNAseq data from GTEx](../../Data_Sources/GTEx_Project/GTEx_RNAseq_Data.md) to try to link genes associated with LDL to specific tissues.

The results are shown below:

![bar_plot](https://github.com/user-attachments/assets/04d531cf-becb-436d-8cbe-24c428f2a55f)


indicating that genes that increase LDL tend to be expressed in the liver. This is consistent with known biology: the liver plays a major role in clearing excess LDL from the bloodstream.

## Reproducing Analysis

To reproduce this analysis, run the [LDL GWAS Analysis Script][mecfs_bio.analysis.ldl_analysis].
