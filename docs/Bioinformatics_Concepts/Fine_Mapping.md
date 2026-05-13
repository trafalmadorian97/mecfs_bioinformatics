

# Fine Mapping

## GWAS and Causal Variants

We typically visualize the results of a GWAS as a Manhattan plot, which shows the profile of univariate-regression p values across the genome. 
The plot below generated from the [DecodeME](../Data_Sources/DecodeME.md)[@genetics2025initial] summary statistics provides an example:

![decode_me_manahattan](https://github.com/user-attachments/assets/149ea85f-c71b-47ea-b208-4dcfc70fa195)

A primary purposes of a GWAS is identification of causal genetic variants.  These causal variants can inform us about the biological processes underlying the phenotype we are studying.

The naive approach is to label the most significant genetic variants  (i.e. the highest points on the Manhattan plot, which are sometimes called "lead variants") as causal.  This approach omits a key consideration: linkage disequilibrium.

## Linkage Disequilibrium
_Linkage disequilibrium_ (LD) refers to statistical dependence between genetic variants. LD is central to statistical genomics.

Some facts about LD:

- LD decays as the distance between variants increases, because the odds of an intervening [recombination event](https://en.wikipedia.org/wiki/Genetic_recombination) correspondingly increase. However, due to the complex structure of Eukaryotic DNA, the odds of recombination events are non-uniform across a chromosome. Thus, the rate of LD decay with genomic distance is not constant.
- Genetic variants that are relatively recent tend to have low frequency in the population, and thus low LD with all other variants, regardless of distance.


As an illustrative example, here is a plot of the absolute value of the correlation between genetic variants in a region of chromosome 1.  This plot was generated from the [UK Biobank LD matrices stored on AWS OpenData](https://registry.opendata.aws/ukbb-ld/).


![ld_example_plot](https://github.com/user-attachments/assets/a05681d5-91f3-4b89-8023-d3d50a22b8bd)

Consistent with the facts above, we observe irregularly spaced LD blocks. 


## Effect of LD on Significance

LD patterns can obscure the true causal variant. To illustrate the point, Wang and Huang[@wang2022methods] provide the toy example illustrated below:


![gwas_wrong_causal_variant](https://github.com/user-attachments/assets/db71a06f-68b5-4410-bc18-3fb3c34c67c9)

In the example, the left and right variants are causal, while the central variant is not.  The left and right variants are correlated with the central variant, but not one another.  As a result, the raw GWAS signal is misleading: if we were to select the candidate causal variant according to p value, we would select the central variant.

Wang and Huang's example reflects a phenomenon seen in real GWAS.  For instance, in a large scale UK Biobank fine-mapping study of 49 traits, Weisbrod et al.[@weissbrod2020functionally] report that "Only 39% of the 2,225 \[likely-causal\] SNPs were also lead GWAS SNPs".



## Fine Mapping

Wang and Huang[@wang2022methods] point out that fine mapping can be thought of as "an exercise to disentangle the effect of LD from the GWAS data".  The goal is to de-convolve the GWAS signal to find the true causal variants.

### Fine Mapping as Bayesian Regression

Fine mapping is often formulated as a Bayesian linear regression variable selection problem.  In the case of a continuous phenotype, we write

$$
y= \sum_{j\in J} \beta_j x_j + \epsilon
$$

where 

- $y$ is the phenotype,
- $J$ is the set of genetic variants at a locus of interest,
- $x_j$ is an individual's allele for variant $j$,
- $\beta_j$ is the effect of variant $j$ on the phenotype,
- $\epsilon$ accounts for the part of the phenotype not attributable to genetic variants at the locus.

We assume a Bayesian prior over the vector $\beta\in\mathbb{R}^J$, and then estimate a posterior over $\beta$ conditional on GWAS data.


Roughly speaking, fine mapping methods conclude that a SNP is causal if there is a high posterior probability that $\beta_i$ is significantly different from zero.  That is:

$$
P(\lvert\beta_i\rvert\gg 0 |D) \approx 1.
$$

where $D$ is the GWAS data. In fine mapping terminology, a causal SNP is one has a high posterior probability of contributing significantly to the phenotype.


It is natural to ask why, when so many other bioinformatic techniques use frequentist statistics, Bayesian statistics has come to dominate fine-mapping.  There are two main reasons.


- First, Bayesian statistics allows a precise description of complex forms of uncertainty.  Suppose that we observe significant marginal effects from 10 SNPs.  Suppose the 10 SNPs are all in mutual tight LD (i.e. they are highly correlated).  Many frequentist statistical techniques would formulate one null hypothesis per SNP, and would be unable to reject any of them.  Effectively, the frequentist techniques would provide no information.  Bayesian methods, however, could conclude with confidence that at least one of the 10 SNPs is causal, but express uncertainty over which. 
- Second, Bayesian statistical techniques allow the incorporation of external information in the form of a prior.  As a simple example, note that a SNP that changes an amino acid in a protein is much more likely to have a biological effect than a SNP that replaces one codon with a [synonym](https://en.wikipedia.org/wiki/Synonymous_substitution).




### Model miss-specification error in fine-mapping

todo



## Links

- [Primer on Fine Mapping by Ran Cui](https://www.youtube.com/watch?v=pglYf7wocSI)

- [Talk on SUSIE Fine Mapping Method by Matthew Stephens](https://www.youtube.com/watch?v=QL3YawgTPhc)