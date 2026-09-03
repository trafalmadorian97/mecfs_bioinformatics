---
tags:
- SuSiE
hide:
- toc
---
# Polyfun Chr1 173.5M-174.5M


## Methodology

The extend the results from my [earlier](../SUSIE/a_Chr1_173M_174M_Locus.md) SUSIE[@wang2020simple] fine-mapping of the DecodeME GWAS-1 signal[@genetics2025initial], I applied SUSIE again, but this time used a prior derived from functional genomic anotations, instead of a uniform prior.


As a linkage disequilibrium reference, I used a [UK Biobank LD matrix hosted on AWS Open Data](https://registry.opendata.aws/ukbb-ld/).  Because this LD reference uses GRCh37 coordinates, I used GWASLab to liftover the DecodeME GWAS-1 summary statistics to GRCh37.

As before, I ran SUSIE four times

- Once with $L=10$,
- Once with $L=2$,
- Once with $L=1$,
- Once with $L=10$ and strict variant filtering.

As before, in my SUSIE runs, I retained palindromic SNPs whose strand orientation GWASLAB was able to determine from allele frequencies in the Thousand Genomes Project, and discarded other palindromic SNPs.

### Prior

For this analysis, I used the precomputed prior provided by the authors of PolyFun[@weissbrod2020functionally][^prior_note]. The PoltFun authors created this prior by running 


[//]: # (To narrow the [DecodeME]&#40;../../../../../Data_Sources/DecodeME.md&#41;[@genetics2025initial] GWAS-1 signal, I [fine-mapped]&#40;../../../../../Bioinformatics_Concepts/Fine_Mapping.md&#41; the hit on chromosome 1 using [SUSIE]&#40;https://stephenslab.github.io/susieR/&#41;[@wang2020simple].)

[//]: # ()

[//]: # ()
[//]: # ()
[//]: # (I ran SUSIE 4 times:)

[//]: # ()
[//]: # (- Once with $L=10$,)

[//]: # (- Once with $L=2$,)

[//]: # (- Once with $L=1$,)

[//]: # (- Once with $L=10$ and strict variant filtering.)

[//]: # ()
[//]: # ($L$ refers to the maximum number of credible sets that can found by SUSIE.  A lower $L$ corresponds to increased regularization, since it decreases the ability of SUSIE to use extra credible sets to fit noise.  Weissbrod et al.[@weissbrod2020functionally] observe that setting $L$ to 1 protects against mismatch between the LD reference population and the GWAS population, because when $L=1$, SUSIE no longer depends on the LD matrix.   They also observe that when $L=2$, even though SUSIE still depends on the LD matrix, empirically it tends to be robust to moderate levels of population mismatch. I thus used the $L=1$ and $L=2$ runs to evaluate whether population mismatch could be influencing SUSIE's results.)

[//]: # ()
[//]: # ()
[//]: # ("Variant filtering" refers to removal of outlier variants according to a [Kriging]&#40;https://en.wikipedia.org/wiki/Kriging&#41;-based likelihood ratio test.  Zou et al.[@zou2022fine] propose this filtering strategy to mitigate instability in SUSIE due to mismatch between the LD and GWAS populations.  In the first three runs above, I filter variants with a likelihood ratio &#40;$\mathrm{LR}$&#41; and absolute $z$ score greater than 2, [consistent with the SUSIE documentation]&#40;https://stephenslab.github.io/susieR/reference/kriging_rss.html&#41;.  In the final run I instead filter variants with $\mathrm{LR}\ge 2$ and $|z|\ge 1$, to evaluate the sensitivity of the results to the filtering threshold.)

[//]: # ()
[//]: # (In my SUSIE runs, I retained palindromic SNPs whose strand orientation GWASLAB was able to determine from allele frequencies in the Thousand Genomes Project, and discarded other palindromic SNPs.)

[//]: # ()
[//]: # ()
[//]: # (## Results)

[//]: # ()
[//]: # (In all 4 runs, SUSIE found a single diffuse credible set.  Moreover, this credible set contained the same 86 variants in all four runs, as illustrated in the UpSet plot below:)

[//]: # ()
[//]: # ()
[//]: # ({{ png_embed&#40;"docs/_figs/decode_mechr1_173500000_174500000_palindromes_keep_upset_plot.png", alt="upset_chrom_1"&#41; }})

[//]: # ()
[//]: # ()
[//]: # (The next figure illustrates the SUSIE results for $L=10$. It is representative.)

[//]: # ()
[//]: # ({{ png_embed&#40;"docs/_figs/decode_mechr1_173500000_174500000_palindromes_keep_susie_stackplot.png", alt="chr1_stackplot"&#41; }})

[//]: # ()
[//]: # (- The top panel is a heatmap in which pixel $&#40;i,j&#41;$ is colored according to the squared correlation between variants $i$ and $j$.  The heatmap reveals the local linkage disequilibrium &#40;LD&#41; structure in the vicinity of the GWAS hit, which is a determinant of SUSIE's results when $L>1$.)

[//]: # ()
[//]: # (- The second panel shows a local Manhattan plot.)

[//]: # ()
[//]: # (- The third panel shows the SUSIE posterior inclusion probability &#40;PIP&#41;.)

[//]: # ()
[//]: # (- The bottom panel shows genes in the region of the GWAS hit. )

[//]: # ()
[//]: # (Overall, SUSIE has returned a diffuse signal in a region with a number of plausible genes.  This makes it unclear which genes deserve follow-up investigation.)

[//]: # ()
[//]: # ()
[//]: # (The expandable table below lists the full SUSIE results for the $L=10$ case.)

[//]: # ()
[//]: # ({{ markdown_table&#40;"docs/_figs/decode_mechr1_173500000_174500000_palindromes_keep_susie_base_convert_cs_to_markdown.mdx", title="Variant List"&#41; }})

[//]: # ()


[^prior_note]: That is, I used the first approach [listed on the PolyFun wiki](https://github.com/omerwe/polyfun/wiki/1.-Computing-prior-causal-probabilities-with-PolyFun).