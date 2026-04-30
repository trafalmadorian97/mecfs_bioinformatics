# Curated Gene Set Analysis



I aimed to run a [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md)-based[@de2015magma] competitive gene set analysis of DecodeME, but suspected that DecodeME would have insufficient statistical power for a fully-unbiased hypothesis-free study.  Therefor, I curated a list of about 50 gene sets from [MSigDB](https://www.gsea-msigdb.org/gsea/index.jsp) selected according to current leading theories of the pathogenesis of ME/CFS.  The list of curated gene sets can be found  {{ api_link("here", "mecfs_bio.constants.curated_gene_set_collections.curated_mecfs_gene_sets") }}. 

The results of MAGMA analysis using these curated gene sets are plotted below:



{{ plotly_embed("docs/_figs/decode_me_gsa_curated_gene_sets_full_magma_bar_plot/magma_gene_set_plot.html", id="DecodeME_Curated_Geneset")
}}



Disappointingly, none of the gene sets meet the Bonferroni-corrected significance threshold. It is interesting to observe the T-cell receptor signalling gene set is the closetest to being significant.  This is consistent with T-cell-centric theories of ME/CFS pathogenesis[@edwardsproposed].


A reasonable next step might be to boost statistical power by pooling multiple studies.  Methods like GenomicSEM[@grotzinger2019genomic] could be used to combine DecodeME both with biobank GWAS of ME/CFS and with distinct but genetically correlated conditions (e.g. multisite pain[@johnston2019genome]).