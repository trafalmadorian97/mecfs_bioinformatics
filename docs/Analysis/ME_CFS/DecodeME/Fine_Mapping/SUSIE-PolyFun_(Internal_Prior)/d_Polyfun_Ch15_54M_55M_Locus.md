---
tags:
  - SuSiE
---
# Chr15 54M-55M

I applied internal-prior PolyFun[@weissbrod2020functionally] [SUSIE](../../../../../Bioinformatics_Concepts/SUSIE.md)[@wang2020simple] [fine-mapping](../../../../../Bioinformatics_Concepts/Fine_Mapping.md) to the [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] Chromosome 15 GWAS-1 signal, using the same methodology as was applied to the [chromosome 1 locus](a_Polyfun_Chr1_173M_174M_Locus.md).



## Results

### Comparison of configurations

The first UpSetPlot below shows that all 4 SUSIE configurations found the same set of variants.


{{
static_img_embed("docs/_figs/decode_me_polyfun_explain_l2_sldsc_priorchr15_54500000_55500000_palindromes_keep_polyfun_upset_all_cs_variants.png",
alt="upset plot for chrom 15")
}}


The second UpSetPlot compares the minimal set of variants needed to achieve 50% PIP across the 4 configurations.  Again, this set of variants is identical across all 4 SUSIE configurations.


{{
static_img_embed("docs/_figs/decode_me_polyfun_explain_l2_sldsc_priorchr15_54500000_55500000_palindromes_keep_polyfun_upset_cs50_variants.png",
alt=" 50 percentile upset plot for chrom 15")
}}

### Detailed Fine mapping results

At this locus, all three priors (uniform, external, internal) produce very similar results.  The plot and table below provide  details.



{{
susie_polyfun_internal_external_explain_plot("docs/_figs/decode_me_polyfun_prior_comparison_chr15_54_925_638_l10_prior_comparison_plot_svg.svg")
}}




{{
susie_polyfun_internal_external_data_table(src="docs/_figs/decode_me_polyfun_prior_comparison_chr15_54_925_638_l10_prior_comparison_table.parquet",
id="chr15_polyfun_susie_table")
}}


These results indicate that at the chromosome 15 locus, there are no relevant high-impact annotations in the S-LDSC baseline model that could tilt either the external or internal priors away from the uniform prior. With all three priors, the most likely variant is **15:55158922:A:G**.