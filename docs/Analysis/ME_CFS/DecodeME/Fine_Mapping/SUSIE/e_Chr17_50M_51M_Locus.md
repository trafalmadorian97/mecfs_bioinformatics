---
hide:
- toc
---
# Chr17 50M-51M


## Methodology

To narrow down the [DecodeME](../../../../../Data_Sources/DecodeME.md)[@genetics2025initial] GWAS-1 signal, I [fine-mapped](../../../../../Bioinformatics_Concepts/Fine_Mapping.md) the GWAS-1 hit on chromosome 17 using SUSIE[@wang2020simple]. The parameters of my runs were the same as for the [chromosome 1 locus](a_Chr1_173M_174M_Locus.md).


## Results

Each of the 4 SUSIE runs produced a single credible set with the same 63 variants, as illustrated by the UpSet plot below:


{{ png_embed("docs/_figs/decode_mechr17_50000000_51000000_palindromes_keep_upset_plot.png", alt="chr17_upset") }}


The distribution of probability mass across these variants was slightly different across the runs, as illustrated by the following UpSet plot, which is restricted to variants of $PIP>0.01$.



{{ png_embed("docs/_figs/decode_mechr17_50000000_51000000_palindromes_keep_upset_plot_pip001.png", alt="chr_17_pip001") }}


The stackplot below, for the $L=10$ case, illustrates the typical PIP pattern across variants.

{{ png_embed("docs/_figs/decode_mechr17_50000000_51000000_palindromes_keep_susie_stackplot.png", alt="chr17_stackplot") }}

The credible set at this locus is focused upstream of the CA10 gene.  This would seem to make a CA10 a good candidate for further investigation.  

The detailed variants results are in the expandable table below:

{{ markdown_table("docs/_figs/decode_mechr17_50000000_51000000_palindromes_keep_susie_base_convert_cs_to_markdown.mdx", title="Variant List (L=10)") }}

