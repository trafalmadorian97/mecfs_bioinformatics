# MAGMA
 

I applied [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to the Neale lab GWAS of self-reported CFS in the UK Biobank.  I used GTEx bulk RNAseq reference data[@gtex2020gtex].  Since this study's phenotype is based on a single-question self-report, it is less precise than the detailed phenotyping of DecodeME.  Nevertheless, it is still interesting and worthwhile to analyze this data.


Results are plotted below:


{{ plotly_embed("docs/_figs/neale_lab_cfs_magma_bar_plot/magma_gene_set_plot.html", id="cfs-neale-magma", caption="Results of applying MAGMA to the Neale Lab GWAS of self-reported CFS in UK Biobank.  y-axis is negative log of p value.  x-axis corresponds to GTEx tissue type.") }}


As might be expected given the low case number and noisy phenotype, there are no significant tissues.  Nevertheless, it is interesting to observe that the tissues with the lowest p-values (brain_cortex, etc) roughly correspond to the significant tissues found when [MAGMA was run on DecodeME](../DecodeME/MAGMA_DecodeME_Analysis.md)

