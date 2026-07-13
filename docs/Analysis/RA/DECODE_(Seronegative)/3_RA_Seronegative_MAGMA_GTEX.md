# MAGMA (GTEx)

I applied [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to the DECODE meta-GWAS of seronegative rheumatoid arthritis[@saevarsdottir2022multiomics]. 


## MAGMA Gene Analysis

I used MAGMA's SNP-wise mean model to run gene-level analysis.


In this step:

- Data from the 1000 genomes projects was downloaded from the [MAGMA website](https://cncr.nl/research/magma/) and used as a linkage disequilibrium reference.
- [Build 151 of dbSNP](https://hgw2.soe.ucsc.edu/cgi-bin/hgTables?hgsid=2912494930_cRufLdpdc1ynRc2sCM3g1WGAWAgH&hgta_doSchemaDb=hg19&hgta_doSchemaTable=snp151Flagged
  ) was used to assign RSIDs to SNPs.
- Magma's default proximity-based rules were used to assign SNPs to genes.


The plot below shows the results of MAGMA gene-level analysis

{{ plotly_embed("docs/_figs/decode_ra_seronegative_magma_gene_manhattan_plot.html", id="magma-ra-neg-gene-manhattan-gene-manhattan-window") }}

It is interesting to note the p-value spike in the HLA region of chromosome 6, consistent with immune pathogenesis. Thus even though there [S-LDSC failed to identify heritability enrichment in immune cells](2_RA_Seronegative_S_LDSC.md), gene-level analysis does point to immune etiology.

Contrasting the above graph with the [corresponding one for seropositive RA](../DECODE_(Seropositive)/3_RA_Seropositive_MAGMA_GTEX.md) we observe both conditions produce spikes in the HLA region of chromosome 6, but the seropositive RA Manhattan plot has more significant genes outside of chromosome 6, and so is less concentrated at that locus.

## MAGMA Gene Property Analysis



I next applied [MAGMA's](../../../Bioinformatics_Concepts/MAGMA_Overview.md) gene property analysis module to seropositive RA.  This step combined the gene analysis results above with tissue-specific RNA expression values from the [GTEx project](../../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex].

The results are plotted below:

{{ plotly_embed("docs/_figs/decode_ra_seronegative_magma_bar_plot/magma_gene_set_plot.html", id="seropositive-ra-magma-gtex-tissue-bar-plot") }}


In this plot, the y-axis corresponds to $-\log_{10}(p)$ values, the x-axis corresponds to tissue type (only tissues with the smallest p values are shown), and bars are colored according to whether their p value meets Bonferroni-corrected significance threshold.

Interestingly, and consistent with the earlier S-LDSC analysis, no tissue types are significant. This may be a consequence of the lack of strong signal outside of the HLA region, in contrast to seropositive RA.

[//]: # (2026-07-11 17:42:40 [debug    ] Figure asset decode_ra_seronegative_magma_gene_manhattan_plot copied to docs/_figs/decode_ra_seronegative_magma_gene_manhattan_plot.html.)
[//]: # (2026-07-11 17:42:40 [debug    ] Figure asset decode_ra_seronegative_ldsc_heritability_markdown copied to docs/_figs/decode_ra_seronegative_ldsc_heritability_markdown.mdx.)
[//]: # (2026-07-11 17:42:40 [debug    ] Figure asset decode_ra_seronegative_hba_magma_plot_extracted copied to docs/_figs/decode_ra_seronegative_hba_magma_plot_extracted.html.)
[//]: # (2026-07-11 17:47:51 [debug    ] Directory figure asset decode_ra_seronegative_magma_bar_plot copied to docs/_figs/decode_ra_seronegative_magma_bar_plot.)