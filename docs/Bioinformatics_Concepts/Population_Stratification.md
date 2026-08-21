# Population Stratification

## Causal privilege

Suppose that we run an epidemiological study with the aim of understanding the effect of variable A on variable B.  We detect a statistical association between A and B. Generally speaking, there are three possible contributors to this association.


- **Forward causality**: A causes B.


``` mermaid
graph LR
A(A) --> B(B);
```


- **Reverse causality**: B causes A.


``` mermaid
graph LR
B(B) --> A(A);
```


- **Confounding**: A and B are both caused by a third factor C.


``` mermaid
graph LR
C(C) --> A(A);
C --> B(B);
```


Determining what part of the association to attribute to each of these possibilities is a challenge that frustrates much of traditional epidemiological research[@hernan2010causal].


A key advantage of genetic studies over traditional epidemiological research is that genetics is causally privileged.  To be more precise, suppose that A is a person's genotype and B is a phenotype of interest observed after birth.  In general:

- A person's genotype is fixed at birth, so reverse causality can be ruled out.
- Most kinds of environmental effects do not affect a person's genotype, so environmental confounding can be ruled out.


Thus a genetic association is much more likely to reflect a forward causal effect than a general epidemiological association.

This causal privilege is a significant advantage, but it does not mean that genetic studies are free of causal inference considerations.  One such consideration is population stratification.


## Types of stratification


### Genetic population stratification

Genetic population stratification occurs when the population under study contains multiple sub-populations, and mating with a sub-population has historically been much more common than mating across sub-populations.  Genetic population stratification can create very-long range correlations between genetic variants (linkage disequilibrium).  Normally, linkage-disequilibrium in humans decays to zero at a distance of a few megabases, and does not cross chromosomal boundaries.  However, genetic population stratification changes this.  For example SNP $p$ on chromosome 1 and SNP $q$ on chromosome 2  may both be more common in a sub-population than in the general population.  Thus having SNP $p$ makes you more likely to be member of the sub-population, which increases your odds of having SNP $q$. Thus SNP $p$ and SNP $q$ are correlated, despite being on different chromosomes. The concept is illustrated by the causal diagram below.



``` mermaid
graph LR
C(Sub-population) --> A(p);
C --> B(q);
```







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
