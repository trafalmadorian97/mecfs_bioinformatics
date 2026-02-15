

# Fine Mapping

## GWAS and Causal Variants

In a GWAS, one runs a single univariate regression per genetic variant.  The results are visualized as a Manhattan plot, which shows the profile of p values across the genome. 
For instance, see this plot generated from the [DecodeME](../Data_Sources/DecodeME.md)[@genetics2025initial] summary statistics:

![decode_me_manahattan](https://github.com/user-attachments/assets/149ea85f-c71b-47ea-b208-4dcfc70fa195)

One of the primary purposes of a GWAS is to identify genetic variants that are causal for the phenotype of interest.  These causal variants can provide insights into the biological processes underlying the phenotype.

It would be natural to take the most significant genetic variants  (i.e. the highest points on the Manhattan plot, which are sometimes called "lead variants") and conclude that they are causal.  This approach omits a key consideration: linkage disequilibrium.

## Linkage Disequilibrium
"Linkage disequilibrium" (LD) refers to statistical dependence between genetic variants. LD is central to nearly every aspect of statistical genomics.

Some facts about LD:

- LD decays as the distance between variants increases, because the odds of an intervening recombination event correspondingly increase. However, due to the complex structure of Eukaryotic DNA, the odds of recombination events are highly non-uniform across a chromosome. Thus, the rate of LD decay with genomic distance is not constant.
- Genetic variants that are relatively recent tend to have low frequency in the population, and thus low LD with all other variants, regardless of distance.


As an illustrative example, here is a plot of the absolute value of the correlation between genetic variants in a region of chromosome 1.  This plot was generated from the [UK Biobank LD matrices stored on AWS OpenData](https://registry.opendata.aws/ukbb-ld/).


![ld_example_plot](https://github.com/user-attachments/assets/a05681d5-91f3-4b89-8023-d3d50a22b8bd)

Consistent with the facts above, we observe irregularly-spaced blocks of high LD.  Moreover, the high LD blocks are not uniform, but contain low-correlation genetic variants, perhaps the result of recent mutation.


## Effect of LD on Significance

LD patterns can obscure the true causal variant. To illustrate the point, Wang and Huang[@wang2022methods] provide the toy example illustrated below:


![gwas_wrong_causal_variant](https://github.com/user-attachments/assets/db71a06f-68b5-4410-bc18-3fb3c34c67c9)

In the example, the left and right variants are causal, while the central variant is not.  The left and right variants are correlated with the central variant, but not one another.  As a result, the GWAS signal is misleading: if we were to select the candidate causal variants according to p value, we would select the central variant.

Wang and Huang's example is not a purely academic construction.  It reflects a phenomenon seen in real GWAS.  For instance, in a large scale UK Biobank fine-mapping study of 49 traits, Weisbrod et al.[@weissbrod2020functionally] report that "Only 39% of the 2,225 \[likely-causal\] SNPs were also lead GWAS SNPs".



## Fine Mapping

Wang and Huang[@wang2022methods] note that fine mapping can be thought of as "an exercise to disentangle the effect of LD from the GWAS data".  The goal is to identify the true variants by reversing the effect of LD.

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

We state a Bayesian prior over the vector $\beta\in\mathbb{R}^J$, and then estimate a posterior over $\beta$ conditional on GWAS data.


Roughly speaking, fine mapping methods will conclude that a SNP is causal if there is a high posterior probability that $\beta_i$ is significantly different from zero.  That is:

$$
P(\lvert\beta_i\rvert\gg 0 |D) \approx 1.
$$

Where $D$ is the GWAS data. In fine mapping terminology, a causal SNP is one has a high posterior probability of contributing significantly to the phenotype.


It is natural to ask why, when so many other bioinformatic techniques use frequentist statistics, Bayesian statistics has come to dominate fine-mapping.  There are two main reasons.


- First Bayesian statistics allows a precise description of complex forms of uncertainty.  Suppose that we observe significant marginal effects from 10 SNPs.  Suppose the 10 SNPs are all mutual tight LD (i.e. they are highly correlated.  Many conservative frequentist statistical techniques would formulate one null hypothesis per SNP, and would be unable to reject any of them.  Effectively, the frequentist techniques would provide no information.  Bayesian methods, however, could conclude with confidence that at least one of the 10 SNPs is causal, but express uncertainty over which. 
- Bayesian statistical techniques allow the incorporation of external information in the form of a prior.  As a simple example, note that a SNP that changes an amino acid in a protein is much more likely to have a biological effect than a SNP that replaces one codon with a [synonym](https://en.wikipedia.org/wiki/Synonymous_substitution).




### Model miss-specification error in fine-mapping

todo



# Links

[Primer on Fine Mapping by Ran Cui](https://www.youtube.com/watch?v=pglYf7wocSI)


