# Heritability
## Definition
Heritability is defined in the context of an additive genetic model.  Let $Y$ denote the phenotype of interest.  Let $G$ denote the genetic component of the phenotype, and let $E$ denote the environmental component.
We assume the contributions of genetic and environmental effects are additively separable, so that:

$$
Y=E+G.
$$

Assume also that $E$ and $G$ are uncorrelated.

Under this assumption, define heritabiity as

$$
h^2:= \frac{\mathrm{Var}(G)}{\mathrm{Var}(Y)}.
$$


Note that in the uncorrelated-additive model above, $h^2$ is equal to the [coefficient of determination](https://en.wikipedia.org/wiki/Coefficient_of_determination) in the linear regression of $Y$ on $G$.



$h^2$ thus expresses the proportion of phenotypic variation attributable to genetics.

## Note on assumption
The assumption that $G$ and $E$ are uncorrelated may be unrealistic in some contexts.  This is problematic for two reasons:


1.  If $\mathrm{Cov}(G,E)\ne 0$, then $h^2$ loses much of its interpretive value.  For instance, it is no longer the coefficient of determination of a regression.
2. According to Visscher et al., if $\mathrm{cov}(G,E)\ne 0$, then  statistical methods for estimating $h^2$ will tend to produce inflated estimates.


As I understand it, practitioners typically argue that any correlation between $G$ and $E$ is not too strong, so that the above effects are not too severe, and the heritability concept remains a useful.







## References

- Balding, David J., Ida Moltke, and John Marioni, eds. Handbook of statistical genomics. John Wiley & Sons, 2019. (Chapter 15)
- Visscher, Peter M., William G. Hill, and Naomi R. Wray. "Heritability in the genomics era—concepts and misconceptions." Nature reviews genetics 9.4 (2008): 255-266.