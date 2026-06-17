
# MAGMA

I  applied gene-level [MAGMA](../../../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to the GWAS-by-subtraction[@demange2021investigating; @huang2024gwas] residual of [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] minus multisite pain[@johnston2019genome].  For my initial analysis, I used upstream and downstream window sizes of 0. As usual, I used the European subset of the thousand genomes project as a linkage disequilibrium reference. The gene-level results are plotted below.

{{ plotly_embed("docs/_figs/decode_me_minus_johnston_ols_magma_gene_manhattan_plot.html", id="subtraction-gene-man", caption="Gene-level Manhattan plot showing the results of standard MAGMA applied to the GWAS-by-subtraction residual of DecodeME minus multisite pain") }}

It is of interest to compare this plot to the [one produced from the original DecodeME summary statistics](../../../../ME_CFS/DecodeME/g_MAGMA_DecodeME_Analysis.md).  Relative to that plot, we observe that the MAGMA signal has been noticeably dampened on Chromosome 12 and toward the end of Chromosome 1.  In contrast, there is a newly significant cluster of genes lead by RFTN2 significant on Chromosome 2.