# MAGMA IBD Analysis
After applying [MAGMA to the DecodeME summary statistics](../DecodeME/MAGMA_DecodeME_Analysis.md), I decided to apply it to summary statistics of a number of other conditions whose biology is better understood.  The aim here was to assess the reliability of MAGMA.  If MAGMA recapitulates known aspects of the pathophysiology of other diseases, this increases our confidence in the results of applying MAGMA to ME/CFS.

## Data
I decided to start with Inflammatory bowel disease (IBD).

I used summary statistics from [Liu et al.'s 2023 cross-ancestry meta-GWAS of inflammatory bowl disease](https://www.nature.com/articles/s41588-023-01384-0). Because my reference data is from the European population, I extracted only the summary statistics corresponding to non-Finish Europeans, excluding those summary statistics from other populations. I passed these summary statistics through the same MAGMA pipeline I applied to DecodeME.

##  Results
As before, I ran the gene analysis set, followed by the gene-set analysis step.  In the gene analysis step, I used [tissue-specific RNAseq data from GTEx](../../Data_Sources/GTEx_project/GTEx_RNAseq_Data.md) to try to link genes associated with IBD to specific tissues.  The results are shown in the bar plot below:

![bar_plot](https://github.com/user-attachments/assets/75bcdd92-5fbd-4722-a2d7-5a019a04fb31)


Of note, the associations with colon and intestinal tissue make sense. These findings increase my confidence in MAGMA.

On the other hand, the association with  lung and whole blood is less obvious. These associations could be spurious, reflecting a limitation of the MAGMA method.  On the other hand, they could be real, albeit surprising, insights into IBD pathology.
