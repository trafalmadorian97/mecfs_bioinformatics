# Colocalization

Colocalization is a biostatistical technique that tests whether two GWAS share a causal variant at a locus.  

Because it evaluates commonality between GWAS signals, colocalization is similar in spirit to [genetic correlation](Genetic_Correlation.md) techniques like [CT-LDSC](Cross_Trait_LDSC.md).  However, whereas genetic correlation evaluates global overlap, colocalization focuses on a single locus.  Two GWAS may strongly colocalize at one locus, but not at another. Colocalization techniques facilitate nuanced understanding of where two traits share biological pathways and where they do not.


Colocalization is often paired with [Mendelian Randomization](Mendelian_Randomization.md), another local technique.


## Bayes factor colocalization

Bayes factor colocalization[@giambartolomei2014bayesian] is simple, but requires assumptions.  It forms the basis for more advanced colocalization methods.



### Data generating model

Assume initially that both traits are quantitative.  Assume the following generative model.


- There are $Q$ genetic variants at the locus of interest.
- There are two GWAS.  There are $N_1$ participants in the first GWAS and $N_2$ participants in the second GWAS.
- $y_k\in\mathbb{R}^{N_k}$ is the vector of phenotypes for participants in GWAS $k$.
- $G_k \in \mathbb{R}^{N_k \times Q}$ is the genotype matrix for GWAS $k$ at the locus.
- $g_{k,i}\in\mathbb{R}^{N_k}$ is the $i$th column of $G_k$.
- $\beta_{k}\in\mathbb{R}^{Q}$ is the vector of true (joint) variant effect sizes for trait $k$.
- $\hat\beta_{k,i}\in\mathbb{R}$ is the GWAS marginal regression coefficient of variant $i$ in GWAS $k$.
- $V_{k,i}\in\mathbb{R}$ is the squared standard error of $\hat\beta_{k,i}$.
- $\epsilon_k\in\mathbb{R}^{N_k}$ is the part of $y_k$ not determined by variants at the locus.  $\epsilon_k$ includes both the environment and other loci.
- $\sigma^2_k>0$ is the variance of the components of $\epsilon_k$.


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

The prior probability of category vector $c\in\mathcal{C}^Q$ is then:

$$
\begin{align}
\prod_{i=1}^Q p_{c_i}. \label{P_def}
\end{align}
$$


Conditional on variant $i$ being causal for trait $k$, we assume the following prior on effect  sizes:

$$
\begin{align}
\beta_{k,i} | (i \text{ is causal for }k) \sim \mathcal{N}(0,\omega_k) \label{effect_prior}
\end{align}
$$ 

for some $\omega_k>0$.


### The single causal variant assumption


Let $\mathcal{A}$ denote the event that at most one variant is causal for each of the two traits.  We next define the function $\overline{S}$, whose domain is the set of category vectors corresponding to $\mathcal{A}$.  Let $\overline{S}(c):=(S_1,S_2)$, where $S_1$ is the index of the causal variant for trait 1, or 0 if trait $1$ has no causal variant, and $S_2$ is the index of the causal variant for trait $2$, or 0 if trait $2$ has no causal variant. A pair $(S_1,S_2)$ is called a configuration. There are $(Q+1)^2$ configurations.
 
Bayes factor colocalization operates conditional on  $\mathcal{A}$, the assumption that at most one variant is causal.  This assumption is unrealistic for some loci, but permits significant simplification.


### Hypotheses


We define 5 hypothesis, which partition the configurations:

- $H_0$ corresponds to the configuration $\mathcal{S}_0=\{ (0,0) \}$.  There are no causal variants at the locus.
- $H_1$ corresponds to the configurations $\mathcal{S}_1=\{ (i,0) : i \ge  1 \}$.  There is a causal variant for trait 1 but not trait 2.
- $H_2$ corresponds to the configurations $\mathcal{S}_2=\{ (0,j) : j \ge  1 \}$.  There is a causal variant for trait 2 but not trait 1.
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
&=\frac{P(S)}{p_0^Q} & \text{by (\ref{P_def}),}\\
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

We compute the probability of each of the hypotheses.  We begin with unconditional probability, then move to probability conditional on $\mathcal{A}$.

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
&=\begin{cases}
\frac{1}{Z} & h=0\\
\frac{Q\rho_1}{Z} & h=1\\
\frac{Q\rho_2}{Z} & h=2\\
\frac{Q(Q-1)\rho_1\rho_2}{Z} & h=3\\
\frac{Q\rho_{1,2}}{Z} & h=4\\
\end{cases} \label{cond_h_prob} \\
\end{align}
$$

where we have defined $Z:=1+Q\rho_1 +Q\rho_2 +Q(Q-1)\rho_{1}\rho_2 + Q\rho_{1,2}$ and where in $(\ref{cond_h_prob})$ we divided the numerator and denominator by $p_0^Q$.




### Marginal GWAS regression coefficients

By the principles of ordinary linear regression, the marginal GWAS regression coefficient $\hat\beta_{i}$ and its standard error $V_{i}$ are given by

$$
\begin{align}
\hat\beta_{i}&=\frac{g_{i}^T y}{\lVert g_{i}\rVert^2}\\
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



where $p$ denotes probability density, and $\{S=i\}$ is the event that $i$ is the sole causal variant, while $\{S=0\}$ is the event that there are no causal variants.  Note that both these events are subsets of the event $\mathcal{A}$ defined above.

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
BF(i)&=\sqrt{1-r_i}\exp \left( \frac{r_i z_i^2}{2} \right). \label{bf_formula}
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
&\exp\left( -\frac{1}{2\sigma^2}( \lVert  y -  g_i \beta_i\rVert^2 -  \lVert  y \rVert^2)  \right)\\
&=\exp\left( -\frac{1}{2\sigma^2}( \lVert y \rVert^2 - 2\beta_i y^T  g_i + \beta_i^2\lVert  g_i \rVert^2 -  \lVert  y \rVert^2)  \right)\\
&=\exp\left( -\frac{1}{2\sigma^2}( \beta_i^2\lVert  g_i \rVert^2 - 2\beta_i y^T  g_i )  \right)\\
&=\exp \left( \frac{\beta_i}{\sigma^2} y^T  g_i- \beta_i^2\frac{\lVert g_i \rVert^2}{2\sigma^2} \right)\\
&=  \exp \left( \beta_i \hat \beta_i \sigma^2 \frac{1}{V_i\sigma^2} -\beta_i^2 \frac{1}{2V_i}\right) & \text{by (\ref{v_rewrite}) and (\ref{beta_rewrite})  }
\end{align}
$$


To produce an expression for the Bayes factor, we need to integrate the above over all possible values of $\beta_i$, weighting by the prior on $\beta_i$ from $(\ref{effect_prior})$. We can apply the Gaussian integral from $(\ref{gaussian_integral})$ with

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


### Interpreting the Bayes Factor

The plot below shows the relationship between the p value of a GWAS variant and the Bayes factor given by $(\ref{bf_formula})$ for a few different values of $\hat\beta$.

{{
static_img_embed("docs/_figs/coloc_bf_vs_pvalue_diagram.svg",alt="image illustrating relationship between p values and Bayes factors", width="100%")
}}


Notable points from the plot:

- In ratio of $\log\mathrm{BF}$ to $-\log p$ approaches 1 as sample size grows.
- However,  $\log\mathrm{BF}$ is not monotonic in $-\log p$.  Initially, $BF$ actually declines as $p$ shrinks.  This effect can be understood as follows: when the data provides evidence against the alternative hypothesis, the alternative hypothesis is more likely at a small sample size (when the data may be unreliable) than at a moderate sample size.
- $\mathrm{BF}$ is not monotonic in $\hat\beta$.  This can be understood as an effect of the prior.  Very large values of $\hat\beta$ are downweighted by the prior.

### Posterior main result

Next, we are interested in how observing the data affects the probability of the hypotheses. 

Define the following quantities

$$
\begin{align}
L_1 &:= \sum_{i=1}^Q B_1(i)\\
L_2 &:= \sum_{i=1}^Q B_2(i)\\
L_{1,2} &:= \sum_{i=1}^Q B_1(i) B_2(i)
\end{align}
$$

Let $D:=(y_1,y_2)$ denote the observed phenotypic data.

Our main result follows.


$$
\begin{align}
P(H_h|D,\mathcal{A})&=\frac{w_h}{\sum_{l=0}^4w_l}
\end{align}
$$

where 

$$
\begin{align}
w_0&:=1\\
w_1&:=\rho_1 L_1\\
w_2&:=\rho_2 L_2\\
w_3&:= \rho_1\rho_2(L_1L_2 -L_{1,2})\\
w_4&:= \rho_{1,2}L_{1,2}
\end{align}
$$


### Proof of main result

First, let $E$ denote any event whose intersection with $\mathcal{A}$ has nonzero probability. We have:


$$
\begin{align}
&P(E|D,\mathcal{A})\\%\frac{P(S|\mathcal{A})}{P(S^0|\mathcal{A})}
&=\frac{P(E\cap \mathcal{A}|D)}{P(\mathcal{A}|D)} & \text{Def of cond prob}\\
&=\frac{p(E\cap \mathcal{A}\cap D)/p(D)}{p(\mathcal{A}\cap D)/p(D)     } & \text{Density ratio} \\
&=\frac{p(E\cap \mathcal{A}\cap D)P(E\cap A)/(p(D)P(E\cap A))}{p(\mathcal{A}\cap D)P(\mathcal{A})/(p(D)P(\mathcal{A}))     } & \text{Multiply and divide}\\
&=\frac{p(D| E, \mathcal{A})P(E\cap A)}{p( D|\mathcal{A})P(\mathcal{A})     }\\
&=\frac{p(D| E, \mathcal{A})P(E| A)}{p( D|\mathcal{A})     } \label{e_prob_breakdown}
\end{align}
$$


Let $S$ denote the random variable whose value corresponds to the one of the $(Q+1)^2$ configurations that describes the true causal structure.  Note that $S$ is only defined on the event $\mathcal{A}$.  Now, apply $(\ref{e_prob_breakdown})$ with the event $E$ equal to $S=s$.


$$
\begin{align}
P(S=s|D,\mathcal{A})&=\frac{p(D|S=s) P(S=s|\mathcal{A})}{p(D|\mathcal{A}|)} \label{s_prob}
\end{align}
$$

Now, we break down the denominator in $(\ref{s_prob})$.

$$
\begin{align}
p(D|\mathcal{A})&=\frac{p(D\cap \mathcal{A})}{P(\mathcal{A})}\\
&=\sum_{s'}\frac{p(D\cap \mathcal{A}\cap\{S=s'\} )}{P(\mathcal{A})}\\
&=\sum_{s'}\frac{p(D\cap \mathcal{A}\cap\{S=s'\} )}{P(\mathcal{A})}\frac{P(\mathcal{A} \cap \{S=s' \})}{P(\mathcal{A} \cap S=s')}\\
&=\sum_{s'} p(D| \mathcal{A},\{S=s' \})P(\{S=s' \}|\mathcal{A})\\
&=\sum_{s'} p(D| \mathcal{A})P(S=s' |\mathcal{A})\\
\end{align}
$$

where we have used the fact that the $(Q+1)^2$ configurations partition $\mathcal{A}$.

In the next step, we divide both the numerator and the denominator by $P(S=S^0|\mathcal{A})p(D|S=S^0)$.


This produces,


$$
\begin{align}
&P(\{ S=s \} |D,\mathcal{A})\\
&=\frac{p(D|S=s)P(S=s|\mathcal{A})}{p(D|\mathcal{A})}\\
&=\frac{p(D|S=s)P(S=s|\mathcal{A})}{\sum_{s'} p(D| S=s')P(S=s' |\mathcal{A})}\\
&=\frac{(p(D|S=s)P(S=s|\mathcal{A})/(P(S=S^0|\mathcal{A})p(D|S=S^0)))}{\sum_{s'} p(D| S=S')P(S=s' |\mathcal{A})/(P(S=S^0|\mathcal{A})p(D|S=S^0))}\\
&=\frac{p(D|S=s)P(S=s|\mathcal{A})}{P(S=S^0|\mathcal{A})p(D|S=S^0)}\frac{P(S=S^0|\mathcal{A})p(D|S=S^0) }{\sum_{s'} p(D| S=s')P(S=s' |\mathcal{A})}
\\
&=\frac{1}{K}\frac{p(D|S=s)P(S=s|\mathcal{A})}{P(S=S^0|\mathcal{A})p(D|S=S^0)}\\
&=\frac{1}{K}\frac{P(S=s|\mathcal{A})}{P(S=S^0|\mathcal{A})}\frac{p(D|S=s)}{p(D|S=S^0)} \label{conf_prob}
\end{align}
$$

Where we have defined $1/K:=\frac{P(S=S^0|\mathcal{A})p(D|S=S^0) }{\sum_{s'} p(D| S=s')P(S=s' |\mathcal{A})}$.

Note that in $\ref{unconditional_prob}$ we have already derived an expression for $\frac{P(S=s|\mathcal{A})}{P(S=S^0|\mathcal{A})}$.

It remains to derive an expression for $\frac{p(D|S=s)}{p(D|S=S^0)}$.

Note that by our assumption of the conditional independence of the two GWAS, we have

$$
\begin{align}
p(D|S=s)&=p(D_1|S_1=s_1) p(D_2|S_2=s_2)\\
\end{align}
$$

Thus $(\ref{conf_prob})$ becomes

$$
\begin{align}
&P(\{ S=s \} |D,\mathcal{A})\\
&=\frac{1}{K}\frac{P(S=s|\mathcal{A})}{P(S=S^0|\mathcal{A})}B_1(s_1)B_2(S_2).
\end{align}
$$

We are now in a position to compute the posterior probabilities

Since each hypothesis is a disjoint union of configurations, we have

$$
\begin{align}
&P(H_h|D,\mathcal{A})\\
&=\sum_{s\in\mathcal{S}_h}P(S=s|D,\mathcal{A})\\
&=\sum_{s\in\mathcal{S}_h} \frac{1}{K}\frac{P(S=s|\mathcal{A})}{P(S=S^0|\mathcal{A})}B_1(s_1)B_2(s_2)\\
&=\frac{1}{K} \rho^h  \sum_{s\in\mathcal{S}_h} B_1(s_1)B_2(s_2)\label{posterior_prob_1}
\end{align}
$$

where 



$$
\begin{align}
\rho^h&:=\begin{cases}
1 & h=0\\
\rho_1 & h=1\\
\rho_2 & h=1\\
\rho_1\rho_2 & h=3\\
\rho_{1,2} & h=4
\end{cases}
\end{align}
$$


Now let us consider $(\ref{posterior_prob_1})$ for each of the possible values of $h$.

For $h=0$, we have


$$
\begin{align}
&P(H_0|D,\mathcal{A})\\
&=\frac{1}{K} \sum_{s\in\mathcal{S}_0} B_1(s_1)B_2(s_2)\\
&=\frac{1}{K}
&=\frac{1}{K} w_0.
\end{align}
$$

For $h=1$, we have

$$
\begin{align}
&P(H_1|D,\mathcal{A})\\
&=\frac{1}{K}\rho_1 \sum_{s\in\mathcal{S}_1} B_1(s_1)B_2(s_2)\\
&=\frac{1}{K}\rho_1 \sum_{1\le s_1 \le Q} B_1(s_1)\\
&=\frac{1}{K}\rho_1 L_1\\
&=\frac{1}{K} w_1.
\end{align}
$$

For $h=2$, we have

$$
\begin{align}
&P(H_2|D,\mathcal{A})\\
&=\frac{1}{K}\rho_2 \sum_{s\in\mathcal{S}_2} B_1(s_1)B_2(s_2)\\
&=\frac{1}{K}\rho_2 \sum_{1\le s_2 \le Q} B_2(s_2)\\
&=\frac{1}{K}\rho_2 L_2\\
&=\frac{1}{K} w_2.
\end{align}
$$


For $h=3$, we have


$$
\begin{align}
&P(H_3|D,\mathcal{A})\\
&=\frac{1}{K}\rho_1\rho_2 \sum_{s\in\mathcal{S}_3} B_1(s_1)B_2(s_2)\\
&=\frac{1}{K}\rho_1 \rho_2 \sum_{1\le s_1,s_2 \le Q, s_1\ne s_2} B_1(s_1) B_2(s_2)\\
&=\frac{1}{K}\rho_1 \rho_2 \left( \sum_{1\le s_1,s_2 \le Q} B_1(s_1) B_2(s_2) - \sum_{1\le s \le Q}B_1(s)B_2(s) \right)\\
&=\frac{1}{K}\rho_1\rho_2 (L_1L_2 - L_{1,2})\\
&=\frac{1}{K} w_3.
\end{align}
$$

Finally, for $h=4$, we have



$$
\begin{align}
&P(H_3|D,\mathcal{A})\\
&=\frac{1}{K}\rho_{12} \sum_{s\in\mathcal{S}_4} B_1(s_1)B_2(s_2)\\
&=\frac{1}{K}\rho_{12} \sum_{1\le s \le Q} B_1(s)B_2(s)\\
&=\frac{1}{K}\rho_{12} L_{1,2}\\
&=\frac{1}{K} w_4
\end{align}
$$


Now, let us study $K$.

$$
\begin{align}
K&= \sum_{s'}\frac{p(D|S=s')}{P(D|S=S^0)}\frac{P(S=s'|\mathcal{A})}{P(S=S^0|\mathcal{A})}\\
&= \sum_{s'}B_1(s'_1)B_2(S'_2)\frac{P(S=s'|\mathcal{A})}{P(S=S^0|\mathcal{A})}\\
&= \sum_{s'\in H_0} B_1(s'_1)B_2(S'_2)\frac{P(S=s'|\mathcal{A})}{P(S=S^0|\mathcal{A})} + \sum_{s'\in H_1} B_1(s'_1)B_2(S'_2)\frac{P(S=s'|\mathcal{A})}{P(S=S^0|\mathcal{A})} \\
&+ \sum_{s'\in H_2} B_1(s'_1)B_2(S'_2)\frac{P(S=s'|\mathcal{A})}{P(S=S^0|\mathcal{A})} + \sum_{s'\in H_3} B_1(s'_1)B_2(S'_2)\frac{P(S=s'|\mathcal{A})}{P(S=S^0|\mathcal{A})} \\
&+ \sum_{s'\in H_4} B_1(s'_1)B_2(S'_2)\frac{P(S=s'|\mathcal{A})}{P(S=S^0|\mathcal{A})}\\
&=w_1+w_2+w_3+w_4
\end{align}
$$


Thus we have


$$
\begin{align}
P(H_h|D,\mathcal{A})&=\frac{w_h}{\sum_{l=0}^4w_l}
\end{align}
$$

as desired




### Interpretation of Bayesian Coloc

The following key quantities are helpful in obtaining an intuitive understanding of Bayes factor colocalization.

$$
\begin{align}
\pi_k(i)&:=\frac{B_k(i)}{L_k}\\
a_k&:=\frac{\rho_kL_K}{1+\rho_kL_k}\\
C&:=\sum_{i=1}^Q \pi_1(i) \pi_2(i)\\
\kappa&:=\frac{\rho_{12}}{\rho_1\rho_2}
\end{align}
$$

to be continued.

## SUSIE colocalization

todo


## Links


[Talk on colocalization at International Statistical Genetics Workshop](https://youtu.be/QNUdAFn1Vas?si=dW0v2QJF0kzZIIea)