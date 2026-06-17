# Genetic Correlation


[As mentioned previously](1_Overview_Subtraction.md), GWAS-by-subtraction[@demange2021investigating; @huang2024gwas] orthogonally decomposes a GWAS trait. As a test of my implementation, I measured the genetic correlation between ME/CFS[@genetics2025initial], multisite pain[@johnston2019genome], and their GWAS-by-subtraction residual $R$ (labeled as _"DecodeME_Minus_Pain_OLS"_  ). The results follow:

{{ markdown_table("docs/_figs/decode_me_minus_pain_ols_genetic_corr_ct_ldsc_corr_agg_markdown.mdx", title="GWAS-by-subtraction genetic correlations") }}

_Key columns: p1: Trait 1; p2: trait 2; rg: genetic correlation; se: standard error_



As expected, the genetic correlation between the residual and multisite pain is statistically indistinguishable from zero.

It is notable that the residual remains very highly correlated with the original DecodeME trait.