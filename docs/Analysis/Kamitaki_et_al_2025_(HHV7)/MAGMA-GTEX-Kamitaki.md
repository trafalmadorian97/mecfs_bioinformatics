---
hide:
  - toc
---
# MAGMA GTEx


I applied [MAGMA](../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to Human Herpesvirus 7 DNA GWAS of Kamitaki et al.[@kamitaki2025genes] using bulk RNAseq data from [GTEx](../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex] as a reference. The gene-property analysis results are plotted below.

{{ plotly_embed("docs/_figs/kamitaki_et_al_2025_hhv7_magma_bar_plot/magma_gene_set_plot.html", id="hhv7-gtex-magma", caption="Result of MAGMA GTEx applied to Kamitaki et al's HHV7 DNA GWAS. The x axis corresponds to tissue type, while the y axis measures MAGMA significance level.") }}

Consistent with the role of the immune system of maintaining HHV7 in a dormant state, the top three tissue types in this analysis are  _Ileum Lymphode Aggregate, Whole Blood_, and _Spleen_, all of which are immune-related.  As noted by Finucate et al., lung tissue in the GTEx dataset also probably contains immune cells[@finucane2018heritability], so it makes sense for lung tissue to be significant. The three other tissue types all gut-related.  Their significant is less clear, but may reflect the role of intestinal tissue as an important site of HHV7 latency[@gonelli2001human].

## Reproducing Analysis

To reproduce the above, run the {{ api_link("HHV7 Analysis Script", "mecfs_bio.analysis.hhv7_analysis") }}.
