# S-LDSC

I applied stratified linkage disequilibrium score regression (S-LDSC)[@finucane2018heritability] to summary statistics from DECODE's meta-GWAS of seronegative rheumatoid arthritis (RA)[@saevarsdottir2022multiomics].


## Reference Data Sources

I used the standard reference datasets prepared by the authors of the S-LDSC method.

- [The GTEx Project](../../../Data_Sources/GTEx_RNAseq_Data.md)
- The Franke lab dataset
- [The Roadmap Epigenetic Project](../../../Data_Sources/Roadmap.md) and the ENCODE project.
- The [Corces et al. ATAC-seq dataset of 13 blood cell types](../../../Data_Sources/Corces_et_al.md).
- The [ImmGen](../../../Data_Sources/Immgen_Project.md) Project
- The Cahoy Mouse Central Nervous System Dataset


## Results

Interestingly, and in contrast to the [analogous analysis for seropositive RA](../DECODE_(Seropositive)/2_RA_Seropositive_S_LDSC.md), there we no significant tissue or cell types. It is unclear what this means.  It could reflect either a mechanistic heterogeneity of patients assigned the diagnostic label "seronegative RA", or could indicate that the key cell type underlying the pathomechanism of seronegative RA is not present in any of the reference datasets I used.