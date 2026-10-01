---
tags:
  - SuSiE
hide:
  - toc
---
# Chr20 47M-48.2M


I applied internal-prior PolyFun[@weissbrod2020functionally] [SUSIE](../../../../../Bioinformatics_Concepts/SUSIE.md)[@wang2020simple] [fine-mapping](../../../../../Bioinformatics_Concepts/Fine_Mapping.md) to the [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] GWAS-1 signal on Chromosome 20, using the same methodology I previously applied to the [chromosome 1 locus](a_Polyfun_Chr1_173M_174M_Locus.md).


The UpSetPlots below illustrate respectively

- The overlap across the four SUSIE runs of all variants found in credible sets, and
- The overlap across the four SUSIE runs of the minimal set of variants required to achieve a total PIP of 50%.



{{
static_img_embed("docs/_figs/decode_me_polyfun_explain_l2_sldsc_priorchr20_47000000_48200000_palindromes_keep_polyfun_upset_all_cs_variants.png",
alt="upset plot for chrom 20 internal prior")
}}



{{
static_img_embed("docs/_figs/decode_me_polyfun_explain_l2_sldsc_priorchr20_47000000_48200000_palindromes_keep_polyfun_upset_cs50_variants.png",
alt="50 PIP upset plot for chrom 20 internal prior")
}}

Unlike the previous loci, but consistent with my [uniform-prior runs at the chromosome 20 locus](../SUSIE/f_Chr20_47M_48M_Locus.md), as well as my [external-prior runs at this locus](../SUSIE-PolyFun_(External_Prior)/f_Chr20_47M_48M_Locus.md) here the results produced by SUSIE depend on the chosen configuration: The $L=10$ and $L=2$ runs produce similar results, distinct from the $L=1$ and strict $L=10$ runs.


### Results ($L=10$)

The plot and table below show the results for the $L=10$ run, which is similar to the $L=2$ run


{{
 susie_polyfun_internal_external_explain_plot("docs/_figs/decode_me_polyfun_prior_comparison_chr20_47_653_230_l10_prior_comparison_plot_svg.svg")
}}


{{ susie_polyfun_internal_external_data_table(src="docs/_figs/decode_me_polyfun_prior_comparison_chr20_47_653_230_l10_prior_comparison_table.parquet",
id="chr20_polyfun_internal_susie_table_l10")}}


### Results ($L=1$)



The plot and table results below show the results for the $L=1$ run, which is similar to the strict $L=10$ run.

{{
susie_polyfun_internal_external_explain_plot(
"docs/_figs/decode_me_polyfun_prior_comparison_chr20_47_653_230_l1_prior_comparison_plot_svg.svg"
)
}}


{{
susie_polyfun_internal_external_data_table(
"docs/_figs/decode_me_polyfun_prior_comparison_chr20_47_653_230_l1_prior_comparison_table.parquet", id="chr_20_l1_table"
)
}}

### Analysis


Given the significance differences between the $L=10$ and $L=2$ SUSIE results on the one hand, and the $L=1$ and strict $L=10$ SUSIE results on the other, it is difficult to know which to credit. As was the case in the [uniform-prior](../SUSIE/f_Chr20_47M_48M_Locus.md) and [external-prior](../SUSIE-PolyFun_(External_Prior)/f_Chr20_47M_48M_Locus.md) SUSIE runs, the $L=10$ and $L=2$ PolyFun-prior runs assigns very high confidence to **20:47743125:C:A**  being a causal SNP, while the $L=1$ and strict $L=10$ run do not weight it at all.


Consistent with our findings at other loci, at this locus the results of the internal-prior runs are generally intermediate between the external-prior and uniform-prior runs. However, this general principle does not universally apply.  For instance, in the $L=1$ run, the top-PIP variant (**20:47526493:T:G**) in the internal-prior run is not highly ranked in either the uniform-prior or external prior runs.