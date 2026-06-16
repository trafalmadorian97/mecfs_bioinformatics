# Genetic Correlation


[As mentioned previously](1_Overview_Subtraction.md), the aim of GWAS-by-subtraction is generate an orthogonal decomposition of a GWAS trait. As a test of my implementation of GWAS-by-subtraction, I measured the genetic correlation between ME/CFS, multisite pain, and the GWAS-by-subtraction residual $R$. The results are below:

{{ markdown_table("decode_me_minus_pain_ols_genetic_corr_ct_ldsc_corr_agg_markdown copied to docs/_figs/decode_me_minus_pain_ols_genetic_corr_ct_ldsc_corr_agg_markdown.mdx", title="GWAS-by-subtraction genetic correlations") }}


As expected, the genetic correlation between the residual and multisite pain is statistically indistinguishable from zero.