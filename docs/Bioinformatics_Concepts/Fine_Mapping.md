

# Fine Mapping

## GWAS and Causal Variants

In a GWAS study, one runs a one univariate regression per genetic variant.  The result of this analysis is typically visualized via a Manhattan plot, showing which univariate regressions produced the most significant regression coefficients.  


See, for instance, this plot from [DecodeME](../Data_Sources/DecodeME.md)

![decode_me_manahattan](https://github.com/user-attachments/assets/149ea85f-c71b-47ea-b208-4dcfc70fa195)

One of the primary purposes of GWAS is to identify genetic variants that play a causal role in the disease or trait of interest.  By finding these causal variants, researchers can derive data-driven insights into the biological processes underlying the phenotype.

It would be natural to take the most significant genetic variants  (i.e. the highest points on the Manhattan plot) and conclude that they are causal.  This, approach, however omits a key consideration: linkage disequilibrium

## Linkage Disequilibrium
"Linkage disequilibrium" (LD) refers to the statistical dependence between genetic variants. LD structure is a central fact that colors nearly every aspect of genomics.

Some facts about LD

- Because genetic variants on the same chromosome are only inherited separately when a recombination event occurs between them, nearby genetic variants tend to have very high LD.  In general, LD decays as the distance between variants increases, because the odds of an intervening recombination event correspondingly increases. However, due the complex structure of Eukaryotic DNA, the odds of recombination events are highly non-uniform across a chromosome.  Thus gross LD structure is non-uniform.
- Genetic variants that are relatively recent on an evolutionary time scale tend to have low frequency in the population, and thus low LD with all other variants.


As an illustrative example, here is a plot of the absolute value of the correlation between genetic variants in region of chromosome 1.


![ld_example_plot](https://github.com/user-attachments/assets/3a396f4b-0da7-4cd6-8345-1d9d6e449688)


## Affect of LD on Significance

LD patterns can obscure the true causal variant. To illustrate the point, Wang et al. [@wang2020simple] provide the following example:

todo


# Links

[Primer on Fine Mapping by Ran Cui](https://www.youtube.com/watch?v=pglYf7wocSI)

