# LDSC

I applied [Linkage Disequilibrium Score Regression](../../../Bioinformatics_Concepts/LDSC.md) (LDSC)[@bulik2015ld] to summary statistics of proteomic GWAS from the European discovery cohort of the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md) (UKBB PPP)[@sun2023plasma]. My aim was to look for evidence of stratification and to estimate [heritability](../../../Bioinformatics_Concepts/Heritability.md).
Note that for many proteomic traits, heritability can be split into two components:



1. A highly concentrated component localized to variants in the vicinity of the gene coding for the protein under study.
2. A diffuse polygenic component spread across the entire genome.


LDSC, because of its modeling assumption of uniform polygenicity, will generally only measure the second heritability component.


As is standard for LDSC analysis, I restricted to the UKBB PPP summary statistics to Hapmap3 variants, and excluded the MHC region.  I used the standard thousand genomes linkage disequilibrium scores provided by the authors of LDSC.


The results are below.


{{ data_table("docs/_figs/ppp_heritability_hapmap_3_eur_discovery_table.parquet", id="ukbb-ppp-ldsc-heritability", caption="SNP heritability computed by LDSC on the UKBB PPP European discovery cohort. Columns: oid: Olink assay ID; gene: name of gene/protein under study;  variant-set: all-variants indicates all hapmap3 variants are included, while cis_excluded indicates that 1-MB region around the gene under study is excluded; h2: LDSC heritability estimate; h2_se: jackknife standard error of heritability; intercept: LDSC intercept; mean_chi2: mean chi squared statistic across included variants; n_snps: number of hapmap3 variants included, which can vary as different variants pass quality control for different proteomic gwas; n_bar: number of participants in cohort with valid proteomic data for the protein under study; p: p value for test that heritability is not zero; ratio: LDSC attenuation ratio, computed only when intercept and mean_chi2 are sufficiently greater than 1." )}}

To validate my results, I compared the all-variant heritabilities reported above and the polygenic heritabilities reported in Supplementary Table 19 of the original Sun et al. [@sun2023plasma] paper.  The Spearman correlation between the two heritabilities was $0.969$, and the mean absolute heritability difference was $0.011$. Thus, there is high level of agreement between my results and those of Sun et al.  The remaining discrepancy is likely due to details of how Sun et al. compute polygenic heritability[^foot].

As for stratification, the intercept and attenuation ratios generally indicate a low stratification effect for the proteins with strong polygenic heritability.  For some of the proteins with weaker polygenic heritability, there is evidence of a stratification effect.


It is interesting to note that the all-variant and cis-excluded heritabilities are generally very similar.  In part, this can be explained by the implementation of LDSC we are using, which, consistent with LDSC's polygenic model, filters SNPs whose $\chi^2$ score exceeds $\max (0.001·N, 80) $.  Thus if proteins have large heritability spikes in their cis regions, these will be filtered out automatically even when we include all variants.


[//]: # (This suggests that, consistent with LDSC model, even when it is provided with summary statistics containing highly localized signal spikes, LDSC will mostly ignore these spikes and compute heritability from a diffuse polygenic signal.  This )






[^foot]:Sun et al. use LDSC to estimate polygenic heritability, but exclude a region around all pQTL, including both cis-pQTL and trans-pQTL.


