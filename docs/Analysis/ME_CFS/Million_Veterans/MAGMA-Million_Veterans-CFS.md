# MAGMA
 

I applied [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to the [Million Veterans](../../../Data_Sources/Million_Veterans_Program.md)[@verma2024diversity] GWAS of Chronic Fatigue Syndrome using bulk RNAseq [GTEx](../../../Data_Sources/GTEx_RNAseq_Data.md) data[@gtex2020gtex] as a reference.  Results are plotted below.



{{ plotly_embed("docs/_figs/million_veterans_cfs_magma_bar_plot/magma_gene_set_plot.html", id="cfs-mv-magma", caption="Results of applying MAGMA to the Million Veterans GWAS of CFS.  y-axis is negative log of p value.  x-axis corresponds to GTEx tissue type.") }}


There are no significant tissues.   This is perhaps to be expected, given the relatively low case count and noisy phenotype definition used here.
