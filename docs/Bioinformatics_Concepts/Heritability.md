# Heritability
## Standard Definition
Heritability is typically defined in the context of an additive genetic model.  Let $Y$ denote the phenotype of interest.  Let $G$ denote the genetic component of the phenotype and $E$ the environmental component.
Assume the components are additively separable, so that:

$$
\begin{align}
Y=E+G. \label{model}
\end{align}
$$

Assume also that the random variables $E$ and $G$ are uncorrelated.

Under this assumption, define heritability as

$$
h^2:= \frac{\mathrm{Var}(G)}{\mathrm{Var}(Y)}.
$$


Note that in the uncorrelated-additive model above, $h^2$ is equal to the [coefficient of determination](https://en.wikipedia.org/wiki/Coefficient_of_determination) in the linear regression of $Y$ on $G$.

To see this, note that if we run a simple linear regression of $Y$ on $G$, then the coefficient of determination is given by the [formula](https://en.wikipedia.org/wiki/Coefficient_of_determination#As_squared_correlation_coefficient):

$$
\begin{align}
\mathrm{cor}(G,Y)^2&=\frac{\mathrm{cov}(Y,G)^2}{\mathrm{Var}(G) \mathrm{Var{Y}}}\\
\mathrm{cor}(G,Y)^2&=\frac{\mathrm{cov}(E+G,G)^2}{\mathrm{Var}(G) \mathrm{Var{Y}}}& \text{ by }(\ref{model})\\
&=\frac{\mathrm{Var}(G)^2}{\mathrm{Var}(G) \mathrm{Var{Y}}} & \text{$G$ and  $E$ uncorrelated}\\
&=\frac{\mathrm{Var}(G)}{\mathrm{Var}(Y)}\\
&=:h^2
\end{align}
$$


$h^2$ thus expresses the proportion of phenotypic variation attributable to genetics.

## Note on assumption
The assumption that $G$ and $E$ are uncorrelated may be unrealistic in some contexts.  This is problematic for two reasons:


1.  If $\mathrm{Cov}(G,E)\ne 0$, then $h^2$ loses much of its interpretive value.  For instance, it is no longer the coefficient of determination of a regression.
2. According to Visscher et al., if $\mathrm{cov}(G,E)\ne 0$, then  statistical methods for estimating $h^2$ will tend to produce inflated estimates.


As I understand it, practitioners typically argue that any correlation between $G$ and $E$ is not strong, so that the above effects are not severe, and the heritability concept remains a useful.
 

## Heritability Depends on the Population Context

The heritability of a trait is highly dependent on the studied population.

For instance, suppose we study a trait in two genetically identical populations $A$ and $B$, but population $B$ lives in a much more variable environment than population $A$, and that this variability affects the trait under study.  Then $\mathrm{Var}(Y_B)>\mathrm{Var}(Y_A)$ .  Consequently, $h^2_A> h^2_B$, despite the genetics of the two populations being the same. 

As a more specific example, note that (a) many recent GWAS have shown that obesity is highly heritable, and that (b) there has been a large increase in obesity rates in recent history, indicating a strong environmental influence.  These two facts do not contradict one another.  The obesity GWASes have studied cross-sectional populations drawn from a single snapshot in time, and have shown that a large portion of the variability in obesity in such cross-sectional populations is attributable to genetics.  On the other hand, if the GWAS samples had included individuals from a range of historical periods, environmental influence would have been proportionately larger, and calculated heritability would have been smaller.

For a popular science book that discusses subtle issues in the definition of heredity, see Harden 2021[@harden2021genetic].

# Alternative definitions of Heritability

## Definition via conditional expectation

If the definition [above](Heritability.md#heritability), we first specified our modeling assumptions, and then defined heritability in the context of those assumptions. It is possible to reverse this order.  One approach is as follows

As before, let the random variable $Y$ denote the phenotypic trait of interest. Let $g$ denote genotype.  Define:

$$
\begin{align}
G&:=\mathbb{E} (Y|g)\\
E&=Y-G\\
h^2&= \frac{\mathbb{Var}(G) }{\mathbb{Var}(Y)}
\end{align}
$$

With this approach, we have

$$
\begin{align}
\mathbb{Cov}(G,E)&=\mathbb{E}(  \mathbb{E}(Y|g) -\mathbb{E}Y   )(  Y- \mathbb{E}(Y|g) )\\
&=0,
\end{align}
$$

where the last line follows from the projection theorem (pg 345 in Grimmet and Stirzaker[@grimmett2020probability]).  Where before we needed to assume $\mathbb{Cov}(E,G)=0$, here this property is automatic.


 - This approach has the **advantage** of its mathematical clarity.  Whereas the standard definition of heritability requires some fairly restrictive assumptions, this definition is applicable to any phenotype representable by a random variable in $L_2$.
 - On the other hand, the conditional expectation definition of heritability has the **disadvantage** of reduced interpretability.  In particular, $G$ now incorporates to everything that can be predicted from the genotype, including artifacts of population stratification.

## Counterfactual heritability





# References

- Balding, David J., Ida Moltke, and John Marioni, eds. Handbook of statistical genomics. John Wiley & Sons, 2019. (Chapter 15)
- Visscher, Peter M., William G. Hill, and Naomi R. Wray. "Heritability in the genomics era—concepts and misconceptions." Nature reviews genetics 9.4 (2008): 255-266.