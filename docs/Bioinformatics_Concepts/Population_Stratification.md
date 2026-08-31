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


A key advantage of genetic association studies is that they are causally privileged.  Let A be a person's genotype and B be a phenotype of interest observed after birth.  In general:

- A person's genotype is fixed at conception, so reverse causality can be ruled out.
- Most kinds of environmental effects do not affect a person's genotype, so environmental confounding can be ruled out.


Thus a genotype-phenotype association is much more likely to reflect a causal effect than a general epidemiological association.

This causal privilege is a significant advantage, but it does not mean that genetic studies are free of the need to account for causal inference considerations.  One such consideration is population stratification[@dattani2022clarifying].


## Types of stratification


### Genetic population stratification

Genetic population stratification occurs when the population under study contains multiple subpopulations, and mating within subpopulations has historically been more common than mating across subpopulations.  Normally, [linkage disequilibrium](Linkage_Disequilibrium.md) in humans decays to zero at a distance of a few megabases, and does not cross chromosomal boundaries.  Genetic population stratification changes this.  For example SNP P on chromosome 1 and SNP Q on chromosome 2 may both be more common in a subpopulation than in the general population due to historical non-random mating.  Thus having  P increases your odds of being a member of the subpopulation, which increases your odds of having Q. P and Q are therefore correlated, despite being on different chromosomes. The concept is illustrated in the causal diagram below.



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

The non-causal association is induced by the backdoor path[^backdoor_note]:

$$
Q \gets \text{Subpopulation} \to P \to \text{Phenotype}.
$$


### Environmental population stratification

It is common for different subpopulations to be exposed to different environments.  These different environments  may differentially affect the phenotype of interest.  This phenomenon is called environmental population stratification.  On its own, environmental population stratification does not confound GWAS.  


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

Even in the extreme case where the phenotype is entirely environmental and does not depend on genetics at all, a combination of genetic and environmental population stratification can induce widespread genotype-phenotype associations.


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



Unfortunately, human population structure is sufficiently complex that it is impossible to fully specify all the subpopulations of which an individual is a member.  Thus, we must use some kind of proxy.  Population genetics research indicates a person's membership in human subpopulations can be well-approximated by the allocation of their genotype to genetic principal components[@price2006principal]. Thus, a strategy to approximately adjust for confounding due to population stratification is as follows[@hoffman2013correcting].


- Let  $X\in\mathbb{R}^{N\times M}$ denote the mean-centered genotype matrix.
- Let $y\in\mathbb{R}^N$ be the mean-centered phenotype vector.
- Let $x_j\in\mathbb{R}^N$ be the $j$th column of $X$.
- Let $\hat\beta_j\in\mathbb{R}$ be the scalar marginal regression coefficient of the $j$th genetic variant.
- Let $\epsilon\in\mathbb{R}^N$ be the random vector of residual environmental and  genetic effects.
- Let $X=USV^T$ be the singular value decomposition of $X$.  Thus $U,V^T \in\mathbb{R}^{N\times N}$ are orthogonal matrices, and $S\in\mathbb{R}^{N\times N}$ is a diagonal matrix.
- Let $q\in\mathbb{Z}_{++}$ be the number of principle components we retain.  Let $U_{1:q}\in\mathbb{R}^{n\times q}$ be the matrix formed from the first $q$ columns of $U$.
- Let $\sigma^2_e>0$ denote the scale of the residual effects.

The PC-controlled marginal GWAS regression for genetic variant $j$ is then described by the following model:

$$
\begin{align}
y &= \hat\beta_j x_j + U_{1:q} \omega + \epsilon\\
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

$\hat\beta_j$ is retained as the marginal GWAS-effect estimate, while $\omega$ is discarded as a nuisance parameter. In this way, we estimate the marginal GWAS association of variant $j$ while controlling for genetic principal components and thus approximately controlling population structure. 



### LMMs

Linear mixed models (LMMs) are another popular method to control for population stratification. Here, I will explain them following Hoffman's derivation[@hoffman2013correcting], which clarifies their connection to the previous controlling-for-PCs approach. Using the same notation as above, consider the following model for the marginal gwas effect of variant $j$.

$$
\begin{align}
y &= \hat\beta x_j + R\gamma + \epsilon \label{gamma_form}\\
\epsilon  &\sim \mathcal{N}(0, \sigma^2_e I)\\
\gamma & \sim \mathcal{N}(0, \sigma^2_\gamma I ),
\end{align}
$$

where
 

- $R :=U_{1:q}S_{1:q,1:q} \in \mathbb{R}^{N \times q}$  
- $\gamma \in \mathbb{R}^q$ is the represents the population structure effect.  Note that in contrast to $\omega$ above, here we put a Bayesian prior on $\gamma$.


We fit this model via maximum likelihood, integrating out $\gamma$:

$$
\begin{align}
\hat\beta_j &= \operatorname*{argmax}_{\hat\beta_j} \int \mathcal{N}(y | \hat\beta_j x_j + R\gamma, \sigma^2_e I ) \mathcal{N}(\gamma|0, \sigma^2_\gamma)  \, \mathrm{d}\gamma
\end{align}
$$


Comparing this model to the direct PC control model of the previous section, the following points are salient:

- In both cases, we us principal components to control for subpopulation membership.
- In the direct PC-control model, we are limited in the number of principal components we can include.  Including too many may result in a model where the number fit parameters approaches or exceeds $N$, the number of study participants, resulting in non-uniqueness of the solution or poor conditioning.  In the LMM, because of the Bayesian prior on $\gamma$, we face no such restriction.
- Because of the presence of $S$ in the definition of $R$, in the LMM more variable principal components can have a larger effect on the phenotype.  In contrast, in the direct PC-control model all components are treated equally.


While the formulation $(\ref{gamma_form})$ is useful for revealing the connection between LMMs and direct PC control, it is not how LMMs are typically writen.  To convert $(\ref{gamma_form})$ to standard LMM form, pick $q=M$ and define $\alpha:= R \gamma$.  By the properties the multivariate normal distribution[^mvnormal_note], 

$$
\begin{align}
\alpha &\sim \mathcal{N}(0, \sigma^2_\gamma RR^T)\\
&= \mathcal{N}(0, \sigma^2_\gamma U S^2 U^T )\\
&= \mathcal{N}(0, \sigma^2_\gamma U SV^T V S U^T )\\\
&=\mathcal{N}(0, \sigma^2_\gamma X X^T)\\\
\end{align}
$$

Thus the restated LMM becomes:

$$
\begin{align}
y &= \hat\beta_j x_j + \alpha + \epsilon\\
\alpha &\sim \mathcal{N}(0,\sigma^2_\gamma K)\\
\epsilon &\sim\mathcal{N}(0, \sigma^2_e I)
\end{align}
$$

Where $K:=XX^T\in\mathbb{R}^{n\times n}$ is called the "genetic relatedness matrix" whose (i,j) entry is a measure of genetic similarity between study participants $i$ and $j$.  This is the standard form that is usually used in presentations of LMMs. 



### LMM Proximal Contamination

While LMMs are effective at controlling for population stratification, if care is not taken they can unduly reduce GWAS statistical power.  This reduction in statistical power can occur for two separate reasons: proximal contamination and ascertainment bias.

We will start by discussing proximal contamination. See the diagram below



``` mermaid
graph LR
C(Subpopulation) --> A(SNP P);
C --> B(G₋ₚ);
A --> D(Phenotype)


classDef normal fill:transparent,stroke:transparent;
classDef conditioned fill:transparent,stroke:#444,stroke-width:2px;
class A,B,C,D normal;
```

Here, P is the causal SNP, while G₋ₚ represents the genome excluding P.   As discussed above, if we use plain regression, population stratification will create false associations between SNPs are the rest of the genome and the phenotype due to confounding by subpopulation.

The LMM approach, especially when $q$ is large, can be understood as controlling for the whole genome as a proxy for controlling for subpopulation.



### LMM Ascertainment Bias


### REGENIE 

todo

### LDSC


todo


[^backdoor_note]: See _Chapter 7: Confounding_ in Hernan and Robins[@hernan2010causal] for a discussion of backdoor paths.

[^mvnormal_note]: See _Section 4.9: Multivariate normal distribution_ from Grimmet and Stirzaker[@grimmett2020probability].

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
