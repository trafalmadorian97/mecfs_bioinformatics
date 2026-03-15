# Univariate MiXeR


Using [software provided by the authors](https://github.com/precimed/mixer), I ran univariate [MiXeR](../../Bioinformatics_Concepts/Mixer.md)[@holland2020beyond] on the DecodeME summary statistics.  

The MiXeR authors recommend running MiXeR 20 times using 20 different random subsets of their standard reference panel of genetic variants.  Here I report my results after 12 runs, since MiXeR runs are quite slow.  I will update this document when the remaining runs finish.

First, here are the key MiXeR output parameters:

| Metric           |          Value |
|:-----------------|---------------:|
| pi (mean)        |    0.00191312  |
| pi (std)         |    0.000173958 |
| sig2_beta (mean) |    4.13238e-05 |
| sig2_beta (std)  |    3.2857e-06  |
| sig2_zero (mean) |    0.966305    |
| sig2_zero (std)  |    0.00339404  |
| h2 (mean)        |    0.162873    |
| h2 (std)         |    0.00432129  |
| nc@p9 (mean)     | 6101.18        |
| nc@p9 (std)      |  554.772       |
| AIC              |   26.5712      |
| BIC              |   36.1355      |


Using Table 2 from the original MiXeR pape[@holland2020beyond] as a reference, we observe that $\pi = 0.0019$ marks ME/CFS as a highly polygenic trait.

Next, consider the power plot generated from the MiXeR model:

![decode_me_mixer_power](../../_figs/decode_me_univariate_mixer_power_plot.png)

MiXeR predicts that an effective sample size of one million would be required to explain slightly more than 20 percent of ME/CFS's heritability with genome-wide significant SNPs.  This conclusion is consistent with ME/CFS being a highly polygenic trait: Since ME/CFS risk is conferred by a large number of weak genetic effects, large statistical power is required to make these weak genetic effects statistically significant.



To assess the goodness of fit of the MiXeR model, we can use a Q-Q plot.

![decode_me_qq](../../_figs/decode_me_univariate_mixer_qq_plot.png)


/// caption
Q-Q plot of DecodeME MiXeR model.  The x axis corresponds to the negative log of the nominal $p$ values, while the y axis corresponds to the negative log the proportion of variants with p values beyond that threshold.  The dotted line is what would be expected under the null hypothesis.  The blue line corresponds to what is actually observed in DecodeME, while the orange line is the predictions of the MiXeR model.
///

This plot suggests that the fitted MiXeR model is a good match for the DecodeME summary statistics.