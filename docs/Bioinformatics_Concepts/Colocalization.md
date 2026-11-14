# Colocalization

Colocalization is a biostatistical technique to evaluate the evidence that two GWAS share a causal variant at a locus.  

Because it evaluates commonality between GWAS signals, colocalization is somewhat similar in spirit to [genetic correlation](Genetic_Correlation.md) techniques like [CT-LDSC](Cross_Trait_LDSC.md).  However, whereas genetic correlation evaluates global overlap, colocalization focuses on a single locus.  Two GWAS may strongly colocalize at one locus, but not at another. Colocalization techniques facilitate nuanced understanding of where two traits share biological pathways and where they do not.


Colocalization is often paired with [Mendelian Randomization](Mendelian_Randomization.md), another local technique.


## Bayes factor colocalization

Bayes factor colocalization[@giambartolomei2014bayesian] is simple, but requires assumptions.  It is used as a basis for more advanced colocalization methods.



### Data generating model

Assume initially that both traits are quantitative.  Assume the following generative model.


- There are $Q$ genetic variants at the locus of interest.
- We have two GWAS.  There are $N_1$ participants in the first GWAS and $N_2$ participants in the second GWAS.
- $y_k\in\mathbb{R}^{N_k}$ is the vector of phenotypes for participants in GWAS $k$.
- $G_k \in \mathbb{R}^{N_k \times Q}$ is the genotype matrix for GWAS $k$ at the locus.
- $g_{k,i}\in\mathbb{R}^{N_k}$ is the $i$th column of $G_k$.
- $\beta_{k}\in\mathbb{R}^{Q}$ is the vector of true variant effect sizes for $k$.
- $\hat\beta_{k,i}\in\mathbb{R}$ i the GWAS marginal regression coefficient of variant $i$ in GWAS $k$.
- $V_{k,i}\in\mathbb{R}$ is the squared standard error of $\hat\beta_{k,i}$.
- $\epsilon_k\in\mathbb{R}^{N_k}$ is the part of $y_k$ not determined by variants at the locus.  $\epsilon_k$ includes both the environment and other loci.
- $c_k$ is the number of covariates in the trait $k$ GWAS model.


Our data generating model is:

$$
\begin{align}
y_k &= G_k \beta_k + \epsilon_k  \label{dg_model} \\
\epsilon_k &\sim N(0,\sigma^2_k I).
\end{align}
$$





We require that $\epsilon_1$ and $\epsilon_2$ are independent.  This requirement is usually satisfied by assuming that there are no overlapping participants between the two studies.



### Category vectors

Variant $i$ is causal for trait $k$ if $\beta_{k,i}\ne 0$.  Use $\mathcal{C}=\{0,1,2,12\}$  to express the 4 disjoint possibilities for any given variant:

- $0$: The variant is not causal.
- $1$: The variant is causal for trait 1 only.
- $1$: The variant is causal for trait 2 only.
- $12$: The variant is causal for both trait 1 and trait 2.

Call $\mathcal{C}^Q$ the set of category vectors.  Each category vector assigns a causal status to each variant.  Assume the causal status of the $Q$ variants are independent and the probabilities of the four possibilities are $p_0,p_1,p_2,p_{1,2}$.

Thus, the prior probability of category vector $c\in\mathcal{C}^Q$ is :

$$
\begin{align}
P(c) = \prod_{i=1}^Q p_{c_i}. \label{P_def}
\end{align}
$$


Conditional on variant $i$ being causal for trait $k$, we assume the following prior on effect  sizes:

$$
\begin{align}
$\beta_{k,i}$ | (i \text{ is causal for }$k$) \sim \mathcal{N}(0,\omega_k) \label{effect_prior}
\end{align}
$$ 

for $\omega_k>0$.


### The single causal variant assumption


Let $\mathcal{A}$ denote the event that at most one variant is causal for each of the two traits.  Define the function $\overline{S}$, whose domain is the set of category vectors corresponding to $\mathcal{A}$.  Define $\overline{S}(c):=(S_1,S_2)$, where $S_1$ is the index of the causal variant for trait 1, or 0 if there is no causal variant, and $S_2$ is the index of the causal variant for trait $2$, or 0 if there is no causal variant. A pair $(S_1,S_2)$ is called a configuration. There are $(Q+1)^2$ configurations.
 
Bayes factor colocalization operates conditional on  $\mathcal{A}$, the assumption that at most one variant is causal.  This assumption is unrealistic for some loci, but permits significant simplification.


### Hypotheses


We define 5 hypothesis, which partition the configurations:

- $H_0$ corresponds to the configurations $\mathcal{S}_0=\{ (0,0) \}$.  There are no causal variants at the locus.
- $H_1$ corresponds to the configurations $\mathcal{S}_1=\{ (i,0) : i \ge  1 \}$.  There is a causal variant for trait 1 but not trait 2.
- $H_2$ corresponds to the configurations $\mathcal{S}_2=\{ (0,j) : j \ge  1 \}$.  There is a causal variant for trait 2 but not trait 1 at the locus.
- $H_3$ corresponds to the configurations $\mathcal{S}_3=\{ (i,j) : i,j \ge  1 ,i\ne j\}$. Both traits have causal variants at the locus, but they are different.
- $H_4$ corresponds to the configurations $\mathcal{S}_4=\{ (i,i) : i \ge  1\}$. The same variant is causal for both traits.

The illustration below shows how the hypotheses partition the configurations.



{{
static_img_embed("docs/_figs/coloc_configuration_grid.svg", width="100%")
}}

### Configuration odds

Define the following probability ratios:

$$
\begin{align}
\rho_1&:=\frac{p_1}{p_0}\\
\rho_2&:= \frac{p_2}{p_0}\\
\rho_{12}&:=\frac{p_{12}}{p_0}.
\end{align}
$$

Let $h\in \{0,1,2,3,4\}$. Let $S\in\mathcal{S}_h$. By $(\ref{P_def})$, we have,


$$
\begin{align}
P(S)&=
\begin{cases}
p_0^Q & h=0\\
p_1 p_{0}^{Q-1} & h=1\\
p_2 p_{0}^{Q-1} & h=2\\
p_1 p_2 p_{0}^{Q-2} & h=3\\
p_{12} p_{0}^{Q-1} & h=4\\
\end{cases} \label{unconditional_prob}
\end{align}
$$



Let $S^0:=(0,0)$. We have


$$
\begin{align}
&\frac{P(S|\mathcal{A})}{P(S^0|\mathcal{A})}\\
&=\frac{P(S \cap \mathcal{A}) / P(\mathcal{A})}{P(S^0 \cap \mathcal{A})/ P(\mathcal{A})}\\
&=\frac{P(S)}{P(S^0)}\\
&=\frac{P(S)}{p_0^Q} & \text{by \ref{P_def},}\\
&=\begin{cases}
1 & h=0\\
\rho_1 & h=1\\
\rho_2 & h=2\\
\rho_1 \rho_2 & h=3\\
\rho_{12} & h=4
\end{cases} & \text{by (\ref{unconditional_prob}).} 
\end{align}
$$


### Prior on the hypotheses

We compute the probability of each of the hypotheses.  We begin with the unconditional probability, then move the probability conditional on $\mathcal{A}$.

There is one configuration corresponding to $H_0$, $Q$ configurations corresponding to $H_1$, $Q$ configurations corresponding to $H_2$, $Q(Q-1)$  configurations corresponding to $H_3$, and $Q$ configurations corresponding to $H_4$.  Thus, by $(\ref{unconditional_prob})$,

$$
\begin{align}
P(H_h)&=\begin{cases}
p_0^Q & h=0\\
Qp_1 p_{0}^{Q-1} & h=1\\
Qp_2 p_{0}^{Q-1} & h=2\\
Q(Q-1)p_1 p_2 p_{0}^{Q-2} & h=3\\
Qp_{12} p_{0}^{Q-1} & h=4\\
\end{cases} &\text{By (\ref{unconditional_prob}).} \label{hy_prob}
\end{align}
$$


$\mathcal{A}$ is the union of these disjoint hypotheses. Thus,

$$
\begin{align}
P(\mathcal{A}) &= P(H_0) + P(H_1) + P(H_2) +P(H_3) +P(H_4)\\
&= p_0^Q + Qp_1 p_{0}^{Q-1} + Qp_2 p_{0}^{Q-1} + Q(Q-1)p_1 p_2 p_{0}^{Q-2} + Qp_{12} p_{0}^{Q-1}.
\end{align}
$$


Now, compute the conditional probability of the hypotheses.


$$
\begin{align}
&P(H_h|\mathcal{A})\\
&=\frac{P(H_h \cap \mathcal{A})}{P(\mathcal{A})}\\
&=\frac{P(H_h)}{P(\mathcal{A})}\\
&=\begin{cases}
\frac{p_0^Q}{P(\mathcal{A})} & h=0\\
\frac{Qp_1 p_{0}^{Q-1}}{P(\mathcal{A})} & h=1\\
\frac{Qp_2 p_{0}^{Q-1}}{P(\mathcal{A})} & h=2\\
\frac{Q(Q-1)p_1 p_2 p_{0}^{Q-2}}{P(\mathcal{A})} & h=3\\
\frac{Qp_{12} p_{0}^{Q-1}}{P(\mathcal{A})} & h=4\\
\end{cases} &\text{By (\ref{hy_prob})}\\
&=\frac{n_h \rho_h}{Z}, \label{cond_h_prob}
\end{align}
$$

where we have defined $\rho_0:=1$ and $Z:=1+Q\rho_1 +Q\rho_2 +Q(Q-1)\rho_{1}\rho_2 + Q\rho_{1,2}$ and 

$$
\begin{align}
n_h:=\begin{cases}
1 & h=0\\
Q &h \in \{1,2,4\}\\
Q(Q-1) & h=3
\end{cases}
\end{align}
$$

and where in $(\ref{cond_h_prob})$ we divided the numerator and denominator by $p_0^Q$.


[//]: # (### Elimination of covariates)

[//]: # ()
[//]: # (The next step is to simplify the data-generating model $&#40;\ref{dg_model}&#41;$ by eliminating the covariate matrix $X$.  The derivations of this section apply separately to the two GWAS, so we suppress the  per-GWAS $k$-subscript.  Let $U\in\mathbb{R}^{&#40;N-c&#41;\times N}$ be a matrix whose rows form an orthonormal basis for the orthogonal complement of the column space of $X$. Define)

[//]: # ()
[//]: # ()
[//]: # ($$)

[//]: # (\begin{align})

[//]: # (\tilde{y}&:=Uy\in\mathbb{R}^{N-c}\\)

[//]: # (\tilde{G}&:=UG\in\mathbb{R}^{&#40;N-c&#41;\times Q}.)

[//]: # (\end{align})

[//]: # ($$)

[//]: # ()
[//]: # (Note that )

[//]: # ()
[//]: # ($$)

[//]: # (\begin{align})

[//]: # (UU^T &= I \in \mathbb{R}^{&#40;N-c&#41;\times &#40;N_c&#41;}\\)

[//]: # (U^TU &= M_X:= I-X&#40;X^TX&#41;^{-1}X^T \in \mathbb{R}^{N \times N}\\)

[//]: # (U\epsilon & \sim \mathcal{N}&#40;0,\omega^2UU^T&#41;=\mathcal{N}&#40;0,\sigma^2I&#41;)

[//]: # (\end{align})

[//]: # ($$)

[//]: # ()
[//]: # (Where $M_X$ is the projection onto the orthogonal complement of the column space of $X$.)

[//]: # ()
[//]: # (Let multiplying equation $&#40;\ref{dg_model}&#41;$ by $U$ produces)

[//]: # ()
[//]: # ($$)

[//]: # (\begin{align})

[//]: # (\tilde{y} &= \tilde{G} \beta + \delta\\)

[//]: # (\delta &\sim \mathcal{N}&#40;0, \sigma^2 I &#41;.)

[//]: # (\end{align})

[//]: # ($$)

[//]: # ()
[//]: # ()
[//]: # (We next relate the )



### Marginal GWAS regression coefficients

By the principles of ordinary linear regression, note that the marginal GWAS regression coefficient $\hat\beta_{i}$ and its standard error $V_{i}$ are given by

$$
\begin{align}
\hat\beta_{i}&=\frac{g_{k,i}^T y}{\lVert g_{i}\rVert^2}\\
V_{i}&= \frac{\sigma^2}{\lVert g_{i}\rVert^2 } \label{v_rewrite}.
\end{align}
$$

where we have suppressed $k$, the study-number subscript.


Combining the above implies that 


$$
\begin{align}
\hat\beta_i&=\frac{ g_i^T  y}{\lVert g_i \rVert^2}\\
&=\frac{ g_i ^T  y}{\sigma^2} V_i \label{beta_rewrite}
\end{align}
$$

### Gaussian integral

The following Gaussian integral will be useful below

$$
\begin{align}
&\int_\mathbb{R} \exp \left(a\beta - b \frac{\beta^2}{2} \right)\mathcal{N}(\beta;0,\omega^2)d\beta\\
&=(1+b\omega^2)^{-\frac{1}{2}} \exp \left( \frac{a^2\omega^2}{2(1+b\omega^2)} \right) \label{gaussian_integral}
\end{align}
$$

It can be proved by standard integration techniques.


### Bayes Factor

As the name suggests, a central quantity in the study of Bayes factor colocalization is the Bayes factor.  Suppressing the study index $k$ for now, the Bayes factor for variant $i$ is defined as

$$
\begin{align}
BF(i):=\frac{p(y|S=i)}{p(y|S=0)}
\end{align}
$$
where $p$ denotes probability density, and $\{S=i\}$ is the event that $i$ is the sole causal variant, while $\{s=0\}$ is the event that there are no causal variants.  Note that both these events are subsets of the event $\mathcal{A}$ defined above.

Define the quantities

$$
\begin{align}
r_i&:=\frac{\omega^2}{\omega^2+V_i}\\
z_i&:=\frac{\hat\beta_i}{\sqrt{V_i}}.
\end{align}
$$

The following result is key:

$$
\begin{align}
BF(i)&=\sqrt{1-r_i}\exp \left( \frac{r_i z_i^2}{2} \right).
\end{align}
$$


**Proof of result**


To prove this, first note that by $(\ref{dg_model})$ we have:

$$
\begin{align}
y | \{S=0,\beta_i \} &\sim \mathcal{N}(g_i \beta_i, \sigma^2 I)\\
y | \{S=0 \} &\sim \mathcal{N}(0,\sigma^2 I).
\end{align}
$$


Now take the ratio of the densities of the two distributions.  We get


$$
\begin{align}
&\exp\left( -\frac{1}{2\sigma^2}( \lVert  y -  g_i \beta_i\rVert^2 -  \lVert  y_i \rVert^2)  \right)\\
&=\exp\left( -\frac{1}{2\sigma^2}( \lVert y \rVert^2 - 2\beta_i y^T  g_i + \beta_i^2\lVert  g_i \rVert^2 -  \lVert  y_i \rVert^2)  \right)\\
&=\exp\left( -\frac{1}{2\sigma^2}( \beta_i^2\lVert  g_i \rVert^2 - 2\beta_i y^T  g_i )  \right)\\
&=\exp \left( \frac{\beta_i}{\sigma^2} y^T  g_i- \beta_i^2\frac{\lVert g_i \rVert^2}{2\sigma^2} \right)\\
&=  \exp \left( \beta_i \hat \beta_i \sigma^2 \frac{1}{V_i\sigma^2} -\beta_i^2 \frac{1}{2V_i}\right) & \text{by (\ref{v_rewrite}) and (\ref{beta_rewrite})  }
\end{align}
$$


To produce an expression for the Bayes factor, we need to integrate the above over all possible values of $\beta_i$, weighting by the prior on $\beta_i$ from $\(ref{effect_prior})$. We can apply the Gaussian integral from $(\ref{gaussian_integral})$ with

$$
\begin{align}
a&=\frac{\hat\beta_i}{V_i}\\
b&=\frac{1}{V_i}
\end{align}
$$

Thus, we have

$$
\begin{align}
BF(i)&=(1+\frac{\omega^2}{V_i})^{-\frac{1}{2}}\exp \left( \frac{\hat\beta_i^2\omega^2 /V_i^2  }{2(1+\omega^2/V_i)}  \right)
\end{align}
$$

Note that

$$
\begin{align}
&(1+\frac{\omega^2}{V_i})^{-\frac{1}{2}}\\
&=(V_i + \omega^2)^{-\frac{1}{2}} V_i^{\frac{1}{2}}\\
&=\sqrt{\frac{V_i}{V_i+\omega^2}}\\
&=\sqrt{1-r_i}
\end{align}
$$

and

$$
\begin{align}
&\frac{\hat\beta_i^2\omega^2 /V_i^2  }{2(1+\omega^2/V_i)} \\
&=\frac{\hat\beta_i^2\omega^2 /V_i  }{2(1+\omega^2)}\\
&=r_iz_i^2/2.
\end{align}
$$


Combining the above, 


$$
\begin{align}
BF(i)&=\sqrt{1-r_i}\exp \left(r_iz_i^2/2 \right)
\end{align}
$$

as desired


**end of proof**

To be continued


## SUSIE colocalization

todo


## Links


[Talk on colocalization at International Statistical Genetics Workshop](https://youtu.be/QNUdAFn1Vas?si=dW0v2QJF0kzZIIea)