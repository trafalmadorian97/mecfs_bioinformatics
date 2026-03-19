# Univariate MiXeR

Using [software provided by the authors](https://github.com/precimed/mixer), I ran univariate [MiXeR](../../Bioinformatics_Concepts/Mixer.md)[@holland2020beyond] on the Johnston et al multsite pain summary statistics[@johnston2019genome].

As recommended by the MiXeR authors, I ran MiXeR 20 times using 20 different subsets of the reference panel of genetic variants.  This serves as a form of [bootstrapping](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)).


The table below lists the key MiXeR output parameters:



--8<-- "docs/_figs/johnston_et_al_pain_univariate_mixer_results_table_as_markdown.mdx"


A polygenicity score of $\pi \approx 0.004$ indicates a very polygenic trait.  A discoverability score of $\sigma^2_\beta \approx 8.85\times 10^{-6}$ indicates very small individual genetic effects.  Compared to ME/CFS, multisite pain appears to be determined by more genetic variants, each of which has a smaller genetic effect.

