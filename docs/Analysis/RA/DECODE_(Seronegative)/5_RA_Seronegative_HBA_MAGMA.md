# MAGMA (HBA)
I applied [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@sey2020computational] to summary statistics from DECODE's meta-GWAS  of seronegative rheumatoid arthritis[@saevarsdottir2022multiomics], using scRNAseq data from the Human Brain Atlas (HBA) as a reference[@siletti2023transcriptomic]. As in my other MAGMA analyses, I sourced linkage disequilibrium reference data from the European subset of the 1000-genomes project. I used a MAGMA gene/cell specificity matrix (i.e. the $E$ matrix in my [notes on MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md) ) prepared as described in Duncan et al.[@duncan2025mapping]


## Results
The results are plotted below


{{ plotly_embed("docs/_figs/decode_ra_seronegative_hba_magma_plot_extracted.html", id="seronegative-ra-hba-magma") }}


There are no significant cell types.  This finding is consistent with the earlier findings from [S-LDSC](2_RA_Seronegative_S_LDSC.md) and [GTEx MAGMA](3_RA_Seronegative_MAGMA_GTEX.md).