# MAGMA EDU Analysis
As another follow-up to the application of [MAGMA to the DecodeME summary statistics](../DecodeME/MAGMA_DecodeME_Analysis.md), I applied MAGMA to the educational attainment study of Lee et al. (2018)[@lee2018gene], to see if the results recapitulated known biology.

## Data
I used summary statistics from [Lee et al.'s GWAS of Educational attainment](https://www.nature.com/articles/s41588-018-0147-3).


##  Results
As before, I ran the gene analysis set, followed by the gene property analysis step.  In the gene property analysis step, I used [tissue-specific RNAseq data from GTEx](../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex] to try to link genes associated with IBD to specific tissues.  The results are shown in the bar plot below:

![bar_plot_edu](https://github.com/user-attachments/assets/d8c4cd10-8b3c-43fa-8246-9b541cac2902)


The significant tissues identified by MAGMA lie exclusively in the brain, consistent with leading theories of the biology of the genetic component of educational attainment.


# How to Reproduce this

To reproduce this analysis run the [Educational Attainment GWAS Analysis Script][mecfs_bio.analysis.lee_educational_attainment_analysis].