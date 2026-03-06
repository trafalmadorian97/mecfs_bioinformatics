---
hide:
- toc
---
# MAGMA GTEx


I applied [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to the lupus GWAS of Bentham et al.[@bentham2015genetic]  using bulk RNAseq data from GTEx[@gtex2020gtex] as a reference.  The gene-property analysis results are plotted below.

{{ plotly_embed("../../../_figs/bentham_2015_lupus_magma_bar_plot/magma_gene_set_plot.html", id="lupus-gtex-magma", caption="Result of MAGMA GTEx applied to Bentham et al.'s lupus GWAS. The x axis corresponds to tissue type, while the y axis measures MAGMA significance level.") }}


The tissue types _small intestine lymphode aggregate, spleen, EBV-transformed_lymphocytes, and whole blood_ are immunological, as befits lupus, an immune-mediated disease.  The relevance of other tissue types, which are gut-related, is less obvious