# MAGMA EDU Analysis
As another follow-up to the application [MAGMA to the DecodeME summary statistics](../DecodeME/MAGMA_DecodeME_Analysis.md), I decided to apply MAGMA to the educational attainment study off Lee Et Al (2018), to see of it recapitulates known biology.

## Data
I used summary statistics from [Lee et al.'s GWAS of Educational attainment](https://www.nature.com/articles/s41588-018-0147-3).


##  Results
As before, I ran the gene analysis set, followed by the gene-set analysis step.  In the gene analysis step, I used [tissue-specific RNAseq data from GTEx](../../Data_Sources/GTEx_project/GTEx_RNAseq_Data.md) to try to link genes associated with IBD to specific tissues.  The results are shown in the bar plot below:

![bar_plot_edu](https://github.com/user-attachments/assets/d8c4cd10-8b3c-43fa-8246-9b541cac2902)


The significant tissues identified by MAGMA lie exclusively in the brain, consistent with leading theories of the biology of the genetic component of educational attainment.


# Reproducing

To reproduce this analysis run the [Educational Attainment GWAS Analysis Script][mecfs_bio.analysis.lee_educational_attainment_analysis]