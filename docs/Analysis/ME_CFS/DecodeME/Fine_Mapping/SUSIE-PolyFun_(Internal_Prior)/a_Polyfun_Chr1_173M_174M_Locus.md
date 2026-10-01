---
tags:
- SuSiE
hide:
- toc
---
# Chr1 173.5M-174.5M


## Methodology

To extend my [earlier](../SUSIE-PolyFun_(External_Prior)/a_Polyfun_Chr1_173M_174M_Locus.md) PolyFun[@weissbrod2020functionally] [SUSIE](../../../../../Bioinformatics_Concepts/SUSIE.md)[@wang2020simple] [fine-mapping](../../../../../Bioinformatics_Concepts/Fine_Mapping.md) of the [DecodeME](../../../../../Data_Sources/DecodeME.md) GWAS-1 signal[@genetics2025initial], I applied PolyFun SUSIE again, but this time used an internal [prior](https://en.wikipedia.org/wiki/Prior_probability) derived the DecodeME summary statistics themselves, rather than an external prior derived from the [UK Biobank](../../../../../Data_Sources/UKBB.md).

In the context of PolyFun SUSIE, the internal and external prior each have their advantages.

- The external prior leverages the statistical power of large UK Biobank GWAS, and is thus less noisy than the internal prior.  Moreover, because certain classes of genetic variants (like evolutionarily conserved variants) carry high heritability across many traits, it is likely that the external prior will contain information highly relevant to any trait under study.
- On the other hand, there may be aspects of the external prior that not applicable to the trait of interest. For instance, it may be that certain annotations are associated with heritability enrichment across most traits, but not in DecodeME.  Using an internal prior mitigates this danger.

I experimented with using an internal prior in my DecodeME fine-mapping to explore this tradeoff.

As a linkage disequilibrium reference, I used a [UK Biobank LD matrix hosted on AWS Open Data](https://registry.opendata.aws/ukbb-ld/).  Because this LD reference uses GRCh37 coordinates, I used [GWASLab](https://github.com/Cloufield/gwaslab) to liftover the DecodeME GWAS-1 summary statistics to GRCh37.

As before, to assess robustness and sensitivity to configuration, I ran SUSIE four times

- Once with $L=10$,
- Once with $L=2$,
- Once with $L=1$,
- Once with $L=10$ and strict variant filtering,

where $L$ refers to the maximum number of credible sets discoverable by SUSIE.

As before, in my SUSIE runs, I retained [palindromic SNPs](../../../../../Bioinformatics_Concepts/Sumstats_Standardization.md#palindromic-variants) and [ambiguous indels](../../../../../Bioinformatics_Concepts/Sumstats_Standardization.md#ambiguous-indels) whose orientation could be determined from allele frequencies in the Thousand Genomes Project.

### Prior Construction

I constructed the internal prior as follows:

1.  I ran $l2$-regularized [stratified linkage disequilibrium score regression](../../../../../Bioinformatics_Concepts/S_LDSC_For_Cell_And_Tissue_ID.md) separately on the odd and even chromosomes of the DecodeME GWAS summary statistics[^annotation_note]. The weight of the $l2$-regularization was determined by within-parity cross validation.  Thus, no information from even chromosomes was used in the odd-chromosome fit, and vice versa.
2. I then used the annotation weights learned from odd-chromosome S-LDSC above to predicts weights for even-chromosome genetic variants based on their annotations, and vice versa.
3. These predicted weights were used as a prior.

The resulting prior upweights genetic variants with functional annotations that are associated with high heritability on opposite parity chromosomes.   The purpose of the odd/even split is to avoid using the same data both to fine map a locus and to learn the prior for that locus.

## Results

### Comparison across runs

I begin by comparing the credible set variants across the $L=1$, $L=2$, $L=10$, and strict $L=10$ runs. The results are plotted in the UpsetPlot below:

{{
static_img_embed("docs/_figs/decode_me_polyfun_explain_l2_sldsc_priorchr1_173500000_174500000_palindromes_keep_polyfun_upset_all_cs_variants.png",
alt="upset plot for chrom 1")
}}

Restricting to the minimal set of variants constituting a total PIP exceeding 50% produces the UpsetPlot:

{{static_img_embed("docs/_figs/decode_me_polyfun_explain_l2_sldsc_priorchr1_173500000_174500000_palindromes_keep_polyfun_upset_cs50_variants.png",
alt="cs50 upset plot for chrom 1")
}}

These results show that at the chromosome 1 locus, SUSIE is insensitive to configuration: we get the same variant set regardless.


### Detailed Fine mapping results


The plot below illustrates the results of $L=10$ SUSIE fine mapping with the uniform prior, the external prior, and the internal prior 


{{
susie_polyfun_internal_external_explain_plot("docs/_figs/decode_me_polyfun_prior_comparison_chr1_174_128_548_l10_prior_comparison_plot_svg.svg")
}}

The table below provides detailed information on $L=10$ SUSIE credible-set variants with and without the PolyFun prior.

{{
susie_polyfun_internal_external_data_table(src="docs/_figs/decode_me_polyfun_prior_comparison_chr1_174_128_548_l10_prior_comparison_table.parquet",
id="chr1_polyfun_susie_table")
}}

Comparing the internal and external prior lift columns in the table above (_lift_int_ and _lift_ext_) we see that at this locus, the two priors boost broadly the same set of variants, but the details differ, resulting in different variants being selected at the top.  Moreover, _lift_ext_ has a higher maximum than _lift_int_, indicating that the external prior is more peaked than the internal prior.  This is consistent with there being a stronger statistical signal in the 15 pooled UK Biobank GWAS used to construct the external prior than in the opposite-parity DecodeME GWAS used to construct the internal prior.







[^annotation_note]: The functional annotations used here come from the baseline model first described in Finucane et al. 2015[@finucane2015partitioning] and extended by Broad Institute researchers.  The version of the baseline model used by PolyFun authors included 187 functional annotations[@weissbrod2020functionally], which cover domains as diverse as evolutionarily conserved regions, QTLs, [epigenetic marks](../../../../../Bioinformatics_Concepts/Epigenetics.md), non-synonymous regions, promoters and enhancers, and more.