---
tags:
- MAGMA
---
# MAGMA HBA

I applied [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to the Pulit et al. meta-GWAS of BMI[@pulit2019meta] using scRNAseq data from the [Human Brain Atlas](../../../Data_Sources/HBA_scRNAseq.md)[@siletti2023transcriptomic] as a reference.


The results are plotted below, first as a static image and then as an interactive plot.

{{ 
static_img_embed("docs/_figs/pulit_et_al_2018_bmi_hba_magma_independent_cluster_plot_svg/hba_magma_fig.svg",alt="static plot of results of HBA MAGMA applied to BMI GWAS", width="100%")   
}}


{{
plotly_embed("docs/_figs/pulit_et_al_2018_bmi_hba_magma_plot_extracted.html",id="hba-magma-bmi")
}}


