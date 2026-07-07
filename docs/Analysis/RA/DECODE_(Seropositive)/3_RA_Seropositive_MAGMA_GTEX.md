# MAGMA-GTEx

I applied [MAGMA](../../../Bioinformatics_Concepts/MAGMA_Overview.md)[@de2015magma] to the DECODE meta-GWAS of seropositive rheumatoid arthritis[@saevarsdottir2022multiomics]. 


## MAGMA Gene Analysis

I used MAGMA's SNP-wise mean model to run gene-level analysis.


In this step:

- Data from the 1000 genomes projects was downloaded from the [MAGMA website](https://cncr.nl/research/magma/) and used as a linkage disequilibrium reference.
- [Build 151 of dbSNP](https://hgw2.soe.ucsc.edu/cgi-bin/hgTables?hgsid=2912494930_cRufLdpdc1ynRc2sCM3g1WGAWAgH&hgta_doSchemaDb=hg19&hgta_doSchemaTable=snp151Flagged
  ) was used to assign RSIDs to SNPs.
- Magma's default proximity-based rules were used to assign SNPs to genes.


The plot below shows the results of MAGMA gene-level analysis

{{ plotly_embed("docs/_figs/decode_ra_seropositive_magma_gene_manhattan_plot.html", id="magma-ra-gene-manhattan-gene-manhattan-window") }}

It is interesting to note the p-value spike in the HLA region of chromosome 6, consistent with immune pathogenesis.

## MAGMA Gene Property Analysis



I next applied [MAGMA's](../../../Bioinformatics_Concepts/MAGMA_Overview.md) gene property analysis module to seropositive RA.  This step combined the gene analysis results above with tissue-specific RNA expression values from the [GTEx project](../../../Data_Sources/GTEx_RNAseq_Data.md)[@gtex2020gtex].

The results are plotted below:

{{ plotly_embed("docs/_figs/decode_ra_seropositive_magma_bar_plot/magma_gene_set_plot.html", id="seropositive-ra-magma-gtex-tissue-bar-plot") }}
In this plot, the y-axis corresponds to $-\log_{10}(p)$ values, the x-axis corresponds to tissue type (only tissues with the smallest p values are shown), and bars are colored according to whether their p value meets Bonferroni-corrected significance threshold,


The tissue types _Whole blood_, _Spleen_, and _EBV-transformed lymphocytes_ make sense, given that RA is an immune-related condition.  The other tissues are less clear.  The presence of lung tissue can be explained by the observation that post-mortem GTEx lung samples may contain blood [@finucane2018heritability].  The gut-related tissue could reflect presence of important immune-related cells in the gut.