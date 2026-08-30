# Population Stratification

## Causal privilege

Suppose we run an epidemiological study to understand the effect of variable A on variable B.  We detect an association between A and B. There are three main possibilities.


- **Causality**: A causes B.


``` mermaid
graph LR
A(A) --> B(B);

classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B normal;
```


- **Reverse causality**: B causes A.


``` mermaid
graph LR
B(B) --> A(A);

classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B normal;
```


- **Confounding**: A and B are both caused by a third factor C.


``` mermaid
graph LR
C(C) --> A(A);
C --> B(B);

classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B,C normal;
```


Determining what part of the association to attribute to each of these possibilities is a challenge that frustrates much traditional epidemiological research[@hernan2010causal].


A key advantage of genetic epidemiology is that it is causally privileged.  Let A be a person's genotype and B be a phenotype of interest observed after birth.  In general:

- A person's genotype is fixed at birth, so reverse causality can be ruled out.
- Most kinds of environmental effects do not affect a person's genotype, so environmental confounding can be ruled out.


Thus a genotype-phenotype association is much more likely to reflect a causal effect than a general epidemiological association.

This causal privilege is a significant advantage, but it does not mean that genetic studies are free of causal inference considerations.  One such consideration is population stratification[@dattani2022clarifying].


## Types of stratification


### Genetic population stratification

Genetic population stratification occurs when the population under study contains multiple subpopulations, and mating within a subpopulation has historically been more common than mating across subpopulations.  Normally, linkage disequilibrium in humans decays to zero at a distance of a few megabases, and does not cross chromosomal boundaries.  Genetic population stratification changes this.  For example SNP P on chromosome 1 and SNP Q on chromosome 2 may both be more common in a subpopulation than in the general population due to historical non-random mating.  Thus having  P increases your odds of being a member of the subpopulation, which increases your odds of having Q. P and Q are therefore correlated, despite being on different chromosomes. The concept is illustrated in the causal diagram below.



``` mermaid
graph LR
C(Subpopulation) --> A(SNP P);
C --> B(SNP Q);


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B,C normal;
```



Suppose now that P has a true causal effect on the phenotype of interest but Q does not.  The long-range correlation between P and Q will produce a GWAS association of Q with the phenotype, creating the false impression of causal GWAS hit in the vicinity of Q.  See below.


``` mermaid
graph LR
C(Subpopulation) --> A(SNP P);
C --> B(SNP Q);
A --> D(Phenotype)


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B,C,D normal;
```

The non-causal association is induced by the backdoor path:

$$
Q \gets \text{Subpopulation} \to P \to \text{Phenotype}.
$$


### Environmental population stratification

It is common for different subpopulations to be exposed to different environments.  These different environments  may differentially affect the phenotype of interest.  This phenomena is called environmental population stratification.  On its own, environmental population stratification does not confound GWAS.  


However, if both genetic and environmental population stratification are present,  environmental stratification can combine with genetic stratification to induce non-causal GWAS associations.  Having both genetic and environmental stratification is common: genetically distinct people often inhabit distinct environments. A possible instance of this kind of confounding is illustrated below


``` mermaid
graph LR
C(Subpopulation) --> D(Environment);
C --> E(SNP)
D --> B(Phenotype)


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class B,C,D,E normal;
```

In the scenario illustrated by the diagram, a non-causal association will be induced between the SNP and the phenotype due to the backdoor path: 


$$
\text{SNP}\gets\text{Subpopulation}\to \text{Environment} \to \text{Phenotype}.
$$

Note that even in the extreme case where the phenotype is entirely environmental and does not depend on genetics at all, a combination of genetic and environmental population stratification can induce widespread genotype-phenotype associations.


## Adjusting for stratification

We have established that although associations in genetic studies are causally privileged, they can still be confounded by population stratification.  So what can be done?  There are a variety of techniques to mitigate the effects of population stratification.


### Controlling for PCs

The classical causal inference strategy to remove confounding is to adjust for the confounder, which in our case is subpopulation membership.  This strategy is illustrated in the causal diagram below, where (as is traditional in the causal inference literature) we draw a box around conditioned variables.


``` mermaid
graph LR
C(Subpopulation) --> D(Environment);
C --> E(SNP)
D --> B(Phenotype)


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class B,D,E normal;
class C conditioned;
```

By conditioning on the subpopulation, we break the non-causal association between the phenotype and the SNP.



Unfortunately, human population structure is sufficiently complex that it is difficult to see how we could ever gather information to fully condition on all possible subpopulation memberships.  Thus, we must use some kind of proxy.  Population genetics research indicates a person's membership in human subpopulations can be well-approximated by the allocation of their genotype to genetic principal components[@price2006principal]. Thus, a strategy to approximately adjust for confounding due to population stratification is as follows


- Let  $X\in\mathbb{R}^{N\times M}$ denote the men-centered genotype matrix.
- Let $y\in\mathbb{R}^N$ be the mean-centered phenotype vector
- Let $x_j\in\mathbb{R}^N$ be the $j$th column of $X$.
- Let $\hat\beta_j\in\mathbb{R}$ be the scalar marginal regression coefficient of the $j$th genetic variant.
- Let $\epsilon\in\mathbb{R}^N$ be the random vector of residual environmental and  genetic effects.
- Let $X=USV^T$ be the singular value decomposition of $X$.  Thus $U.V^T \in\mathbb{R}^{N\times N}$ are orthogonal matrices, and $S\in\mathbb{R}^{N\times N}$ is a diagonal matrix.
- Let $q\in\mathbb{Z}_{++}$ be the number of principle components we retain.  Let $U_{1:q}\in\mathbb{R}^{n\times q}$ be the matrix formed from the first $q$ columns of $U$.
- Let $\sigma^2_e>0$ denote the scale of the residual effects.

The PC-controlled marginal GWAS regression for genetic variant $j$ is then[@hoffman2013correcting] described by the following model:

$$
\begin{align}
y &= \hat\beta x_j + U_{1:q} \omega + \epsilon\\
\epsilon  &\sim \mathcal{N}(0, \sigma^2_e I).
\end{align}
$$

We estimate $\hat\beta_j$ by maximum likelihood:

$$
\begin{align}
\hat\beta_j,\omega &= \operatorname*{argmax}_{\hat\beta_j,\omega} \mathcal{N}(y|\hat\beta x_j + U_{1:q} \omega,\sigma^2_e I)
\end{align}
$$

where $\mathcal{N}(y|\mu,\Sigma)$ denotes the multivariate normal density with mean $\mu$ and covariance $\Sigma$ evaluated at $y$.

$\hat\beta_j$ is retained as the marginal GWAS-effect estimate, while $\omega$ is discarded as a nuisance parameter. In this way, we estimate the marginal GWAS effect while controlling for genetic principal components, which serve as a proxy for population structure. 



### LMMs






### REGENIE 

todo

### Downstream adjustment: LDSC


todo

[//]: # (A key advantage of genetic studies over non-genetic epidemiological studies is that genetic studies are causally privileged.  Specifically, genetic studies benefit from the following advantages:)

[//]: # ()
[//]: # (1.  A person's genes are fixed at birth. Therefore, when we detect an association of a genotype with a phenotype observed later in life, we can be confident there is no reverse causation.  That is, the direction of causality is from the genotype to the phenotype, and not the reverse.)

[//]: # (2. )

[//]: # ()
[//]: # (``` mermaid)

[//]: # (graph LR)

[//]: # (A[Genotype] --> C[Phenotype];)

[//]: # (B[Environment] --> C;)

[//]: # (```)
