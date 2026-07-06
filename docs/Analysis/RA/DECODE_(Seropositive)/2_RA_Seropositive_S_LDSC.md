# S-LDSC

I applied stratified linkage disequilibrium score regression (S-LDSC)[@finucane2018heritability] to summary statistics from DECODE's meta-GWAS of seropositive rheumatoid arthritis[@saevarsdottir2022multiomics].


## Reference Data Sources

I used the standard reference datasets prepared by the authors of the S-LDSC method.

- [The GTEx Project](../../../Data_Sources/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../../Data_Sources/Roadmap.md) and the ENCODE project.
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../../Data_Sources/Corces_et_al.md).
- The [ImmGen](../../../Data_Sources/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset


## Results

### GTEx and Franke lab tissue expression data
The plot and expandable table below show the results of S-LDSC applied to RA using the GTEx and Franke lab gene expression reference datasets. In the plot, the x-axis corresponds to cell type, while the y-axis corresponds to $-\log_{10}(p)$.  Points are colored according to tissue category.  Large points indicate cell/tissue types deemed significant by the Benjamini-Hochberg procedure at an FDR of 0.01[@benjamini1995controlling].  
