---
tags:
  - SuSiE
---
# Chr17 50M-51M

I applied internal-prior PolyFun[@weissbrod2020functionally] [SUSIE](../../../../../Bioinformatics_Concepts/SUSIE.md)[@wang2020simple] [fine-mapping](../../../../../Bioinformatics_Concepts/Fine_Mapping.md) to the [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] GWAS-1 signal on Chromosome 17, using the same methodology I applied to the [chromosome 1 locus](a_Polyfun_Chr1_173M_174M_Locus.md).


### Comparison of configurations

The UpSetPlots below illustrate respectively

- The overlap across the four SUSIE runs of all variants found in credible sets, and
- The overlap across the four SUSIE runs of the minimal set of variants required to achieve a total PIP of 50%.


These plots show that in general, the variant sets selected by the four runs are very similar, though the $L=1$ run differs slightly from the others.


{{
static_img_embed("docs/_figs/decode_me_polyfun_explain_l2_sldsc_priorchr17_50000000_51000000_palindromes_keep_polyfun_upset_all_cs_variants.png",
alt="upset plot for chrom 17")
}}


{{
static_img_embed("docs/_figs/decode_me_polyfun_explain_l2_sldsc_priorchr17_50000000_51000000_palindromes_keep_polyfun_upset_cs50_variants.png",
alt="50 PIP upset plot for chrom 17")
}}


### Detailed Fine mapping results


The plot and table below illustrates the result of the $L=10$ PolyFun-prior SUSIE run.


{{
susie_polyfun_internal_external_explain_plot("docs/_figs/decode_me_polyfun_prior_comparison_chr17_50_237_377_l10_prior_comparison_plot_svg.svg") 
}}


{{
 susie_polyfun_internal_external_data_table(src="docs/_figs/decode_me_polyfun_prior_comparison_chr17_50_237_377_l10_prior_comparison_table.parquet",
id="chr17_internal_polyfun_susie_table")
}}


At this locus the internal prior seems intermediate between the external prior and the uniform prior.  It is not as flat as the uniform prior, but not as peaked as the external prior.  Like the external prior run, the internal prior run favors evolutionarily conserved variants, and assigns the highest PIP to **17:50291040:C:T**.  This is in contrast to the uniform prior run, which assigns a PIP of only about 1 percent to this variant.