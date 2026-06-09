---
hide:
- toc
---
# Chr6 97.5M-99M


## Methodology

To narrow the [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] GWAS-1 signal, I [fine-mapped](../../../../../Bioinformatics_Concepts/Fine_Mapping.md) the second hit on chromosome 6 using SUSIE[@wang2020simple]. The parameters of my SUSIE runs were the same as for the [chromosome 1 locus](a_Chr1_173M_174M_Locus.md).


## Results

The 4 SUSIE runs each produced a single credible set containing 98 variants.  These  credible sets were largely but not entirely overlapping:

{{ png_embed("docs/_figs/decode_mechr6_97500000_99000000_palindromes_keep_upset_plot.png", alt="chr6b_upsetplot") }}


The stacked plot for the $L=10$ run:

{{ png_embed("docs/_figs/decode_mechr6_97500000_99000000_palindromes_keep_susie_stackplot.png", alt="ch6b_stackplot") }}

And for the $L=2$ run:

{{ png_embed("docs/_figs/decode_mechr6_97500000_99000000_palindromes_keep_susie_stackplot_2_credible_set.png", alt="chrb_stackplot_l2") }}

Overall, SUSIE returns a diffuse signal in a region rather far from genes, making the genetic mechanism of action non-obvious. 


The variants in the credibles set for the $L=10$ and $L=2$ runs are listed below:



{{ markdown_table("docs/_figs/decode_mechr6_97500000_99000000_palindromes_keep_susie_base_convert_cs_to_markdown.mdx", title="Variant List (L=10)") }}

{{ markdown_table("docs/_figs/decode_mechr6_97500000_99000000_palindromes_keep_susie_2_convert_cs_to_markdown.mdx", title="Variant List (L=2)") }}

