---
hide:
- toc
---

# MAGMA GTEx

## Tissue Enrichment

[MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] was applied to the UK BioBank GWAS of migraine[@uk2025whole][@GCST90473326_migraine] using bulk RNAseq data from [GTEx](../../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex] as a reference. The gene-property analysis results are plotted below.

{{ plotly_embed("docs/_figs/uk_biobank_2025_migraine_eur_magma_bar_plot/magma_gene_set_plot.html", id="uk-biobank-2025-migraine-eur-magma-bar-plot") }}

The four significant tissues in this analysis were all digestive tissues. This result is consistent with another migraine GWAS which also used the UK BioBank (but also included other individuals), and where the top MAGMA tissues were also all digestive tissues.[@choquet2021new] 

The finding of significant digestive tissues is somewhat surprising, as migraine is generally considered a neurological disease[@ruschel2024migraine]. This MAGMA finding may relate to observations that irritable bowel syndrome is found with increased prevalence in those with migraine.[@wongtrakul2022increased]

A different migraine GWAS in a USA cohort found MAGMA enrichment primarily in brain tissues, but also in the uterus.[@gasperi2025multi]