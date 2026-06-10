---
hide:
- toc
---
# Chr15 54.5M-55.5M


## Methodology

To narrow the [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] GWAS-1 signal, I [fine-mapped](../../../../../Bioinformatics_Concepts/Fine_Mapping.md) the hit on chromosome 15 using SUSIE[@wang2020simple]. The parameters of my SUSIE runs were the same as for the [chromosome 1 locus](a_Chr1_173M_174M_Locus.md).


## Results

Each of the 4 SUSIE runs produced a single credible set containing identical variants, as illustrated by the following UpSet plot:

{{ png_embed("docs/_figs/decode_mechr15_54500000_55500000_palindromes_keep_upset_plot.png", alt="ch15upsetplot") }}


Here is a plot of the SUSIE results for the $L=10$ case.  It is representative of the other runs.


{{ png_embed("docs/_figs/decode_mechr15_54500000_55500000_palindromes_keep_susie_stackplot.png", alt="CH15_susie_results") }}


Compare to the other loci, the SUSIE credible set at this locus is relatively concentrated on a single variant: a G $\to$ A substitution at GRCh37 position 55158922 on chromosome 15.

The detailed results for the $L=10$ run are below.


{{ markdown_table("docs/_figs/decode_mechr15_54500000_55500000_palindromes_keep_susie_base_convert_cs_to_markdown.mdx", title="Variant List (L=10)") }}

