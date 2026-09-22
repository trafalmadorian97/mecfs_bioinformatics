---
tags:
- SuSiE
hide:
- toc
---
# Chr20 47M-48.2M


## Methodology

To narrow the [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] GWAS-1 signal, I [fine-mapped](../../../../../Bioinformatics_Concepts/Fine_Mapping.md) the hit on chromosome 20 using [SUSIE](../../../../../Bioinformatics_Concepts/SUSIE.md)[@wang2020simple]. The parameters of my SUSIE runs were the same as for the [chromosome 1 locus](a_Chr1_173M_174M_Locus.md).


## Results


At the chromosome 20 locus, the $L=10$ run returned three credible sets, while the $L=2$ run returned two.  One of these credible sets is highly concentrated, consisting of a single $\mathrm{PIP}=1.0$ variant, while the others are diffuse.  In contrast, the $L=1$ and strict $L=10$ SUSIE runs return a single diffuse credible set.  It is unclear which of these results is correct.  It would be interesting to re-run SUSIE on this locus using individual-level data.


As expected, the UpSet plot for this locus reveals that the variants found by the $L=1$ and strict $L=10$ runs are identical, while the other runs differ. 



{{ static_img_embed("docs/_figs/decode_mechr20_47000000_48200000_palindromes_keep_upset_plot.png", alt="chr20_upset_plot") }}


### L=10 Result


The Stackplot and variant list for the base $L=10$ run shows the highly concentrated nature of one of the credible sets and the diffuse nature of the others.


{{ static_img_embed("docs/_figs/decode_mechr20_47000000_48200000_palindromes_keep_susie_stackplot.png", alt="chr20_l10") }}



{{ markdown_table("docs/_figs/decode_mechr20_47000000_48200000_palindromes_keep_susie_base_convert_cs_to_markdown.mdx", title="Variant List (L=10)") }}

### Strict L=10 Result

In contrast, the stackplot for the strict $L=10$ run reveals a single credible set


{{ static_img_embed("docs/_figs/decode_mechr20_47000000_48200000_palindromes_keep_susie_stackplot_strict.png", alt="chr20_strict") }}



{{ markdown_table("docs/_figs/decode_mechr20_47000000_48200000_palindromes_keep_susie_strict_convert_cs_to_markdown.mdx", title="Variant List (L=10, Strict)", collapse_threshold=10) }}


This region of chromosome 20 is relatively gene-dense, so that are a number of plausible causal genes.



