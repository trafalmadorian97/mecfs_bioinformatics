# Population Stratification


Suppose that we run an observational study with the aim of understanding the effect of variable A on variable B.  When, in our observational study, we detect a statistical association between A and B, there are typically three possibilities:


- **Forward causality**: A causes B.


``` mermaid
graph LR
A --> B;
```


- **Reverse causality**: B causes A.


``` mermaid
graph LR
B --> A;
```


- **Confounding**: A and B are both caused by a third factor C.


``` mermaid
graph LR
A --> B;
C --> A;
C --> B;
```


Distinguishing between these possibilities is a very tough problem that frustrates much of traditional epidemiological research[@hernan2010causal].


A key advantage of genetic studies over traditional epidemiological research is that genetics is causally privileged.  To be more precise, if
suppose we detect an association between a genotype and a phenotype observed after birth.

-  A person's genotype is fixed at birth, so reverse causality can be ruled out.
- Most kinds of environmental effects do not affect a person's genotype, so environmental confounding can be ruled out.


Thus a genetic association is much more likely to be causal than an association detected in general epidemiology.



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
