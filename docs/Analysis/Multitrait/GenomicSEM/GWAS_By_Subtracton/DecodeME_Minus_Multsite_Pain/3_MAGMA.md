
# MAGMA

I first applied gene-level [MAGMA](../../../../../Bioinformatics_Concepts/MAGMA_Overview.md) to the GWAS-by-subtraction residual of DecodeME minus multisite pain.  For my initial analysis, I used a window size of 0. As usual, I used the European subset of the thousand genomes project as a linkage disequilibrium reference. The results are plotted below.

{{ plotly_embed("assets/base_asset_store/gwas/me_cfs_minus_pain_ols/decode_me_minus_johnston/analysis/plots/decode_me_minus_johnston_ols_magma_gene_manhattan_plot.html", id="subtraction-gene-man", caption="Gene-level Manhattan plot showing the results of standard MAGMA applied to the GWAS-by-subtraction residual of DecodeME minus multisite pain") }}

It is of interest to compare this plot to the [one produced by the original DecodeME summary statistics](../../../../ME_CFS/DecodeME/g_MAGMA_DecodeME_Analysis.md).  Relative to that plot, we observe that the MAGMA signal has been noticeably dampened on Chromosome 12 and toward the end of Chromosome 1.  In contrast, there is a newly significant cluster of genes lead by RFTN2 significant on Chromosome 2.