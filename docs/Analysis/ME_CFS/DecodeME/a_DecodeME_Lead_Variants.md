---
hide:
- toc
---
# Lead Variants

As an initial analysis step, we apply [GWASLab's procedure for extracting lead variants](https://cloufield.github.io/gwaslab/ExtractLead/) to the DecodeME GWAS-1 data.  This procedure groups together significant genetic variants using a sliding-window approach, then reports the most significant variant in each region.

## Table of Variants from GWASLAB

Here is the table of lead variants produced by GWASLab:

{{ markdown_table("docs/_figs/decode_me_gwas_1_lead_variants_markdown.mdx", title="Lead Variants") }}

Note that the genomic coordinates in the table above refer to genome build 38.

## Comparison with lead variants reported in DecodeME
The lead variants reported above by GWASLAB all agree with the variants reported by the DecodeME preprint[@genetics2025initial].  However, in some cases, GWASLAB assigns these lead variants to different genes than the DecodeME preprint.

In particular


-  17:52183006:C:T is assigned to CA10 by both GWASLAB and the DecodeME preprint.
-  20:48914387:T:TA is assigned to both ARFGEF2 and CSE1L by the DecodeME preprint.  GWASLAB assigns only ARFGEF2.
-  The other variants are assigned to different genes by GWASLAB and the DecodeME preprint.  This reflects ambiguity about how to assign GWAS signals to genes.

## Reproducing
To reproduce this analysis, run the {{ api_link("initial DecodeME analysis script", "mecfs_bio.analysis.decode_me_initial_analysis") }}.


