---
hide:
- toc
---
# cis-pQTL MR Analysis


I applied [Mendelian Randomization](../../../Bioinformatics_Concepts/Mendelian_Randomization.md) to DecodeME GWAS-1[@genetics2025initial] using cis-pQTL from the [UK Biobank Pharma Proteomics Project](../../../Data_Sources/UKBB_PPP.md)[@sun2023plasma] as instruments.  I aimed to identify candidate proteins that may be causal in the ME/CFS disease process.


Since each protein in the UKBB PPP has at most one cis-pQTL, I used the Wald Ratio method to estimate the causal effect of the proteins on ME/CFS.

I used the R package `TwoSampleMR`.


**NOTE**: since the controls for the DecodeME project come from the UK Biobank, and since the pQTLs from the UKBB PPP are also derived from analysis of the UK Biobank, we do not actually have two fully independent samples, so we partially violate the assumptions of `TwoSampleMR`.  This could introduce some error. 

## Results

`TwoSampleMR` identified two potentially causal plasma proteins: RABGAP1L, significant at the Bonferroni-corrected level of 0.01, and BTN1A1, which additionally emerges when the threshold is relaxed to 0.05.  Both are listed in the table below and discussed in turn in the following subsections.

{{ markdown_table("docs/_figs/decode_me_gwas_1_mr_cis_mr_markdown.mdx", title="cis-pQTL MR Results") }}

### RABGAP1L

At a significance level of 0.01 (Bonferroni-corrected), `TwoSampleMR` identified one potentially causal plasma protein: RABGAP1L.

The mechanism by which RABGAP1L could affect ME/CFS is not immediately clear.  The[ Candidate Genes supplement to the DecodeME preprint](https://www.pure.ed.ac.uk/ws/portalfiles/portal/533352484/Candidate_Genes.pdf) notes that "RABGAP1L inhibits the entry and/or promotes the expulsion of bacteria or viruses. It is a viral restriction factor. Reduction in RABGAP1L expression would
be expected to enhance susceptibility to bacterial and viral infection, which often precedes
initial ME/CFS symptoms."


To further investigate, I generate two region plots.  The first shows the DECODE ME GWAS signal in region of the RABGAP1L gene:

{{ png_embed("docs/_figs/decode_me_region_plot_rabgap1l_locus/region_plot.png", alt="mecfs_rabgap1l_region_plot") }}


The second shows the same region from the UK Biobank Pharma Proteomics Project GWAS of plasma RABGAP1L levels:

{{ png_embed("docs/_figs/rabgap1l_region_plot_rabgap1l_locus/region_plot.png", alt="rabgap1l_region_plot") }}

From visual inspection, it seems clear that at least the primary signals do not colocalize: the GWAS peaks are clearly in different locations.  However, it could be that there are secondary signals that do colocalize.  A colocalization algorithm may clarify the question.

### BTN1A1

If we relax the significance threshold to 0.05, one additional potentially causal protein is identified: BTN1A1.


This protein appears to be anti-inflammatory, and higher levels are associated with a greater risk of ME/CFS.


Again, to further investigate, I compared region plots in the vicinity of BTN1A1 from DecodeME and from the UK Biobank Pharma Proteomics Project GWAS of plasma BTN1A1 levels.

The first plot is from DecodeME:

{{ png_embed("docs/_figs/decode_me_region_plot_btn1a1_locus/region_plot.png", alt="decode-me-btn1a1-region") }}


The second is from UK Biobank Pharma Proteomics Project:


{{ png_embed("docs/_figs/btn1a1_ukbb_ppp_region_plot/region_plot.png", alt="ukbb_ppp_btn1a1") }}


In this comparison, it seems even more obvious that the primary signals do not colocalize: they are separated by a recombination boundary.


## Follow-up questions

- Does co-localization analysis suggest that the DecodeME causal variant and the cis-pQTL co-localize?

- Do the cis-pQTLs identified above co-localize with trans-pQTLs for other proteins? If so, this could suggest causal protein regulatory chains.
