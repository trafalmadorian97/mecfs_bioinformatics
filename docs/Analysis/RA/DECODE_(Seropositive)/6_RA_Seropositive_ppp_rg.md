# PPP Genetic Correlation

I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate [genetic correlation](../../../Bioinformatics_Concepts/Genetic_Correlation.md) between the DECODE meta-GWAS of seropositive rheumatoid arthritis[@saevarsdottir2022multiomics] and proteomic GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].


Note that because of its assumptions of uniform polygenicity, CT-LDSC mostly measures genetic correlation due to diffuse polygenic effects.  It does not accurately measure genetic correlation due to highly-concentrated locus-specific effects.

## Results

As is standard for LDSC analysis, I restricted to the statistics to Hapmap3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. Because my [previous heritability experiment](../../Proteomics/UKBB_PPP/PPP_LDSC.md) suggested little difference between using all SNPs and excluding the cis-region near the protein of interest, I only ran this experiment with the cis region excluded.

The results are below:

{{ data_table("docs/_figs/seropositive_ra_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-ra-rg", caption="Columns: oid: Olink assay ID; gene: name of gene/protein under study; rg: CT-LDSC genetic correlation estimate; rg_se: jackknife standard error of CT-LDSC genetic correlation; rg_p: p value of test that rg is not zero; gcov: genetic covariance; inter: intercept term in CT-LDSC regression; h2_trait: trait heritability estimate; h2_prot: protein heritability estimate; n_snps: number of hapmap3 variants included; spr: for cases in which multiple rows corresponding to Olink assays of the same protein have been merged into a single row, this gives the maximum spread between the rg values of the merged rows; s_bh: True if the null hypothesis is rejected under the Benjamini-Hochberg procedure at an FDR of 0.05; s_bon: True if the null hypothesis is rejected under the Bonferroni correction at a significance level of 0.05." )}}

Note that the trait heritability reported above differs slightly from [what we found earlier](4_RA_Seropositive_LDSC.md).  This is due to the excluclusion of SNPs whose $\chi^2$ score exceeds $\max (0.001·N, 80) $ in the implementation of LDSC we use for our proteomic genetic correlation analysis.


The proteins with the most significant genetic correlations are consistent with known rheumatoid arthritis biology:

- CCL19 and CCL21 are known to be over-expressed in patients with rheumatoid arthritis [@pickens2011characterization].
- [TNFRSF9](https://en.wikipedia.org/wiki/TNFRSF9) and [PDCD1](https://en.wikipedia.org/wiki/Programmed_cell_death_protein_1) are important immune checkpoint molecules.
- [TNFRSF8](https://en.wikipedia.org/wiki/CD30), also known as CD30, is a receptor expressed by activated lymphocytes.


## Caveats



These results are intriguing, but it is important to keep in mind the limitations of the quantity we are calculating.  A significant genetic correlation tells us that there is a shared genetic architecture between rheumatoid arthritis and the plasma level of a protein of interest. However, genetic correlation does not inform us about the causal relationship between the protein and the trait.

