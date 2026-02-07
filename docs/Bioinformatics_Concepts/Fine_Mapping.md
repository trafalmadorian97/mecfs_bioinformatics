

# Fine Mapping

## GWAS and Causal Variants

In a GWAS, one runs one univariate regression per genetic variant.  The result is visualized via Manhattan plot, showing which univariate regressions produced the most significant regression coefficients.  See, for instance, this plot from [DecodeME](../Data_Sources/DecodeME.md):

![decode_me_manahattan](https://github.com/user-attachments/assets/149ea85f-c71b-47ea-b208-4dcfc70fa195)

One of the primary purposes of GWAS is to identify genetic variants that play a causal role in the phenotype of interest.  By finding these causal variants, researchers can derive data-driven insights into the biological processes underlying the phenotype.

It would be natural to take the most significant genetic variants  (i.e. the highest points on the Manhattan plot) and conclude that they are causal.  This approach omits a key consideration: linkage disequilibrium

## Linkage Disequilibrium
"Linkage disequilibrium" (LD) refers to statistical dependence between genetic variants. LD is a central to nearly every aspect of statistical genomics.

Some facts about LD:

- LD decays as the distance between variants increases, because the odds of an intervening recombination event correspondingly increases. However, due the complex structure of Eukaryotic DNA, the odds of recombination events are highly non-uniform across a chromosome.  Thus the rate of LD decay with genomic distance is not constant.
- Genetic variants that are relatively recent tend to have low frequency in the population, and thus low LD with all other variants, regardless of distance.


As an illustrative example, here is a plot of the absolute value of the correlation between genetic variants in a region of chromosome 1.


![ld_example_plot](https://github.com/user-attachments/assets/a05681d5-91f3-4b89-8023-d3d50a22b8bd)

Consistent with the facts above, we observe irregularly-spaced blocks of high LD.  Moreover, the high LD blocks are not uniform, but contain low-correlation genetic variants, perhaps the result of recent mutation.


## Effect of LD on Significance

LD patterns can obscure the true causal variant. To illustrate the point, Wang and Huang[@wang2022methods] provide a toy example illustrated by the figure below:


![gwas_wrong_causal_variant](https://github.com/user-attachments/assets/db71a06f-68b5-4410-bc18-3fb3c34c67c9)

In the example, the left and right variants are causal, while the central variant is not.  The left and right variants are correlated with the central variant, but not one another.  As a result, the GWAS signal is misleading: if we were to select the candidate causal variant purely on p value, we would select the central variant.



## Fine Mapping

Wang and Huang[@wang2022methods] note that fine mapping can be thought of as "an exercise to disentangle the effect of LD from the GWAS data".  The goal is to identify the true variants by reversing the obscuring effect of LD.

Fine mapping is naturally formulated as a Bayesian linear regression variable selection problem.  In the case of a continuous phenotype, one writes

$$
y= \sum_{j\in J} \beta_j x_j + \epsilon
$$

where 

- $y$ is the phenotype
- $J$ is the set of genetic variants at the locus
- $x_j$ is an individual's allele for variant $j$
- $\beta_j$ is the effect of variant $j$ on the phenotype
- $\epsilon$ accounts for the part of the phenotype not attributable to genetic variants at the locus.



# Links

[Primer on Fine Mapping by Ran Cui](https://www.youtube.com/watch?v=pglYf7wocSI)

