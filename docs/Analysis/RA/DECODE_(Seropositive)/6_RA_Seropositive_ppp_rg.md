# PPP Genetic Correlation

I applied [Cross Trait Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/Cross_Trait_LDSC.md) (CT-LDSC)[@bulik2015atlas] to estimate genetic correlation between the DECODE meta-GWAS of seropositive rheumatoid arthritis[@saevarsdottir2022multiomics] and proteomic GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma].


Note that because of its assumptions of uniform polygenicity, CT-LDSC mostly measures genetic correlation due to diffuse polygenic effects.  It does not accurately measure genetic correlation due to highly-concentrated locus-specific effects.


As is standard for LDSC analysis, I restricted to the statistics to Hapmap3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC. Because my [previous heritability experiment](../../Proteomics/UKBB_PPP/PPP_LDSC.md) suggested little difference between using all SNPs and excluding the cis-region near the protein of interest, I only ran this experiment with the cis region excluded.

The results are below:

{{ data_table("docs/_figs/seropositive_ra_ppp_rg_cis_excluded_display_frame.parquet", id="ukbb-ppp-ldsc-ra-rg", caption="Columns: oid: Olink assay ID; gene: name of gene/protein under study;  variant-set: all-variants indicates all hapmap3 variants are included, while cis_excluded indicates that 1-MB region around the gene under study is excluded; rg: CT-LDSC genetic correlation estimate; rg_se: jackknife standard error of CT-LDSC genetic correlation; rg_p: p value of test that rg is not zero; gcov: genetic covariance; gcov_intercept: intercept term in CT-LDSC regression; h2_trait: trait heritability estimate; h2_protein: protein heritability estimate; n_snps: number of hapmap3 variants included." )}}


The proteins with the most significant genetic correlations are consistent with known rheumatoid arthritis biology:

- CCL19 and CCL21 are known to be over-expressed in patients with rheumatoid arthritis [@pickens2011characterization].
- [TNF inhibitors](https://en.wikipedia.org/wiki/TNF_inhibitor) are a major class of rheumatoid arthritis treatment.
- [TNFRSF9](https://en.wikipedia.org/wiki/TNFRSF9) and [PDCD1](https://en.wikipedia.org/wiki/Programmed_cell_death_protein_1) are important immune checkpoint molecules.


These results are impressive, but it is important to keep in mind the limitations of the quantity we are calculating.  A significant genetic correlation tells us that there is a partially shared genetic architecture between rheumatoid arthritis and the plasma level of a protein of interest. However,

- A significant genetic correlation does not directly inform us about the causal relationship between RA and the protein.
- 