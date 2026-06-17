# MAGMA HBA

I applied Human Brain Atlas[@siletti2023transcriptomic] (HBA) [MAGMA](../../../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] gene property analysis to the GWAS-by-subtraction[@demange2021investigating; @huang2024gwas] residual of [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] minus Johnston et al.'s GWAS of multisite pain[@johnston2019genome].  The results are plotted twice, first as a static plot, then as an interactive plot.



{{ png_embed("docs/_figs/decode_me_minus_johnston_ols_hba_magma_independent_cluster_plot/hba_magma_fig.png", alt="") }}


{{ plotly_embed("docs/_figs/decode_me_minus_johnston_ols_hba_magma_plot_extracted.html", id="subtraction-hba-magma") }}




It is interesting to compare these results to the [HBA MAGMA analysis of the raw DecodeME summary statistics](../../../../ME_CFS/DecodeME/h_MAGMA-HBA-DecodeME.md).  In that analysis, there were three distinct independently significant clusters, one of which was a type of "Eccentric Medium Spiny Neuron" (EMSN). In the present GWAS-by-subtraction analysis only one independently significant cluster remains, and this cluster is an EMSN.  These results suggest that the EMSN supercluster may be important to distinguishing ME/CFS from multisite pain.