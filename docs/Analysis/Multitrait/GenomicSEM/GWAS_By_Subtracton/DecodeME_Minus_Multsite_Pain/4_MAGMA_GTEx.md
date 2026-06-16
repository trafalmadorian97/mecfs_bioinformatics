

# MAGMA GTEx
I next applied [MAGMA's](../../../Bioinformatics_Concepts/MAGMA_Overview.md) gene property analysis module to the GWAS-by-subtraction residual of [DecodeME](../../../Data_Sources/DecodeME.md) minus Johnston et al.'s GWAS of multisite pain.  This step combined the [zero-window gene-level MAGMA results](3_MAGMA.md) with tissue-specific RNA expression values from the [GTEx project](../../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex].  The aim was to identify tissues enriched for genes associated with ME/CFS.  The results are plotted below:


{{ plotly_embed("docs/_figs/decode_me_minus_johnston_ols_magma_bar_plot/magma_gene_set_plot.html", id="subtraction-gene-bar-gtex", caption="Bar plot showing the result of GTEx-based gene property analysis applied to the GWAS-by-subtraction residual of DecodeME minus multisite pain.  The x axis indicates tissue type, while the y axis corresponds to negative log p value.") }}



The results are qualitatively similar to the output of the [MAGMA GTEx analysis of the original DecodeME summary statistics](../../../../ME_CFS/DecodeME/g_MAGMA_DecodeME_Analysis.md).  This is likely because the GTEx tissue types are insufficiently fine-grained to distinguish between the two sets of summary statistics.