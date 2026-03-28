# DecodeME/Multisite Pain 

I applied [Bivariate MiXeR](../../../../Bioinformatics_Concepts/Mixer.md)[@frei2019bivariate] to summary statistics from [DecodeME](../../../../Data_Sources/DecodeME.md)[@genetics2025initial] and Johnston et al.'s GWAS of multisite pain[@johnston2019genome] in an attempt to quantify their joint genetic architecture.

In the [standard workflow suggested by authors](https://github.com/precimed/mixer/blob/499bc4955cfb3a6024311064cf5c9c896879cc83/usecases/mixer_real/MIXER_REAL.job#L25-L31) of bivariate MiXeR, 20 separate bivariate MiXeR models are trained, each on a different random subset of a reference panel of SNPs. This serves as a form of bootstrapping.  Each of these models are then evaluated on the full reference panel.  Unfortunately, MiXeR is memory-hungry, requiring something like 32GB for this evaluation step.  To reduce memory requirements, I evaluated on random subsets of the reference panel instead of the full reference panel.  In the future, I may rent a cloud machine to do full evaluation on the entire reference panel.


The key output parameters from my bivariate MiXeR runs are shown below:


--8<--  "docs/_figs/initial_bivariate_mixer_DecodeME_Multisite_pain_bivariate_mixer_results_table_as_markdown.mdx"

The most important observation from these output parameters is that _best_vs_min_AIC_ and _best_vs_max_AIC_ are both slightly negative.  According to the [README provided by the authors](https://github.com/precimed/mixer/blob/master/README.md#aic-bic-interpretation)  this indicates that there is not strong statistical support for the full bivariate MiXeR model over simpler statistical models. I assume that this finding is a consequence of noise in the DecodeME summary statistics due to the relatively small sample size of the DecodeME study[^aic_training_note].

This weak statistical support for application of bivariate MiXeR means that findings from bivariate MiXeR for this pair of traits need to be taken with a grain of salt.  We can still analyze bivariate MiXeR's out, but this analysis should be regarded as exploratory, not definitive.

With this caveat, the main finding of bivariate MiXeR applied to ME/CFS and Multisite pain is that the two traits share a large polygenic overlap.  This result is illustrated in the standard MiXeR plot below:


![me_pain_bivariate_mixer](https://github.com/user-attachments/assets/145f1a43-6ce2-49d5-8f6f-b64ad13bae55)




[^aic_training_note]: MiXeR computes AIC as part of the model fitting procedure, so AIC is not affected by my decision to evaluate on random subsets of the reference panel instead of the full reference panel.


The conclusion from the Venn diagram is rather striking: as far as bivariate MiXeR can tell, the biological processes that cause ME/CFS are mostly a subset of the biological processes that cause chronic, multisite pain.  Only a relatively small number of genetic variants truly distinguish ME/CFS.



