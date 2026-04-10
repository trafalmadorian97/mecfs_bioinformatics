# Univariate MiXeR

Using [software provided by the authors](https://github.com/precimed/mixer), I ran univariate [MiXeR](../../Bioinformatics_Concepts/Mixer.md)[@holland2020beyond] on the Johnston et al multsite pain summary statistics[@johnston2019genome].

As recommended by the MiXeR authors, I ran MiXeR 20 times using 20 different subsets of the reference panel of genetic variants.  This serves as a form of [bootstrapping](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)).


The table below lists the key MiXeR output parameters:



{{ include_file("docs/_figs/johnston_et_al_pain_univariate_mixer_results_table_as_markdown.mdx") }}


A polygenicity score of $\pi \approx 0.004$ indicates a very polygenic trait.  A discoverability score of $\sigma^2_\beta \approx 8.85\times 10^{-6}$ indicates very small individual genetic effects.  Compared to ME/CFS, multisite pain appears to be determined by more genetic variants, each of which has a smaller genetic effect.

As befits a polygenic trait with weak effects, a very large sample size would be required to explain a significant proportion of SNP heritability.  This is illustrated by the power plot:

![pain_mixer_power](../../_figs/johnston_et_al_pain_univariate_mixer_power_plot.png)


The QQ-plot suggests a good fit between the Mixer model and the multi-site pain GWAS data:

![pain_qq](../../_figs/johnston_et_al_pain_univariate_mixer_qq_plot.png)