// Bayes-factor colocalization (coloc.abf): a self-contained derivation.
// Compile: pixi r typst compile experiments/claude/coloc_abf_derivation/coloc_abf_derivation.typ
// Every identity below is checked numerically by verify_coloc_abf_derivation.py.

#set document(title: "Bayes-factor colocalization: a self-contained derivation")
#set page(paper: "a4", margin: (x: 2.5cm, y: 2.5cm), numbering: "1")
#set text(font: "New Computer Modern", size: 11pt)
#set par(justify: true, leading: 0.62em, spacing: 0.95em)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)", supplement: [Eq.])
#show heading.where(level: 1): set block(above: 1.6em, below: 0.8em)
#show heading.where(level: 2): set block(above: 1.3em, below: 0.7em)
#show link: set text(fill: rgb("#1f4e79"))

#let accent = rgb("#1f4e79")
#let statement-kind = "statement"
#show figure.where(kind: statement-kind): set align(left)
#show figure.where(kind: statement-kind): set block(breakable: true, above: 1em, below: 1em)

#let statement(supplement, title: none, body) = figure(
  kind: statement-kind,
  supplement: supplement,
  numbering: "1",
  outlined: false,
  block(width: 100%, inset: (left: 10pt, y: 3pt), stroke: (left: 1.5pt + accent))[
    #strong[#supplement #context counter(figure.where(kind: statement-kind)).display()]#if title != none [ (#title)]#strong[.]
    #body
  ],
)
#let lemma(title: none, body) = statement([Lemma], title: title, body)
#let theorem(title: none, body) = statement([Theorem], title: title, body)
#let proposition(title: none, body) = statement([Proposition], title: title, body)
#let corollary(title: none, body) = statement([Corollary], title: title, body)
#let remark(title, body) = block(above: 0.9em, below: 0.9em)[_Remark (#title)._ #body]
#let proof(body) = block(above: 0.7em, below: 1em)[_Proof._ #body #h(1fr) $square$]

#let PP = math.op("PP")
#let BF = math.op("BF")
#let ABF = math.op("ABF")
#let bh = $hat(beta)$
#let dist = $tilde.op$
#let gt = $tilde(g)$
#let yt = $tilde(y)$
#let Normal = $cal(N)$
#let Sset = $cal(S)$
#let Aev = $cal(A)$

#align(center)[
  #text(size: 17pt, weight: "bold")[Bayes-factor colocalization] \
  #v(1pt)
  #text(size: 12.5pt)[A self-contained derivation of `coloc.abf`] \
  #v(3pt)
  #text(size: 10pt)[31 October 2026]
]
#v(6pt)

Two genome-wide association studies (GWAS) both show a signal at the same locus. Is that because
one variant influences both traits, or because two different variants, correlated through linkage
disequilibrium (LD), each influence one trait? Colocalization by approximate Bayes factors,
introduced by Giambartolomei et al. [1], answers with a posterior probability over five hypotheses,
using only per-variant summary statistics and no LD matrix.

This note derives the method from first principles. It states the model, derives every formula, and
identifies exactly which steps are exact and which are approximate. Sources are listed at the end,
and no result here depends on reading them. Beyond the derivation, the note gives a
four-number reduction of the posterior (@thm-reduction), which is the clearest way to see how the
method behaves.

#block(width: 100%, fill: rgb("#eef3f8"), inset: 10pt, radius: 3pt)[
  #set math.equation(numbering: none)
  *The algorithm in brief.* Each trait $k in {1,2}$ supplies an effect estimate $bh_(k i)$ and
  standard error $sqrt(V_(k i))$ at each of the same $Q$ variants $i$. Given a prior standard
  deviation $omega_k$ of the causal effect, set
  $
    z_(k i) = bh_(k i) / sqrt(V_(k i)), quad
    r_(k i) = omega_k^2 / (omega_k^2 + V_(k i)), quad
    Lambda_(k i) = 1/2 [ log(1 - r_(k i)) + r_(k i) z_(k i)^2 ],
  $
  $
    L_1 = sum_(i=1)^Q e^(Lambda_(1 i)), quad
    L_2 = sum_(i=1)^Q e^(Lambda_(2 i)), quad
    L_(12) = sum_(i=1)^Q e^(Lambda_(1 i) + Lambda_(2 i)).
  $
  With prior parameters $p_1, p_2, p_(12)$, the posterior probabilities of the hypotheses are
  $PP_h = w_h \/ sum_(ell=0)^4 w_ell$, where
  $
    (w_0, w_1, w_2, w_3, w_4)
    = (thin 1, thin p_1 L_1, thin p_2 L_2, thin p_1 p_2 (L_1 L_2 - L_(12)), thin p_(12) L_(12) thin).
  $
]

= The model <sec-model>

== Data

Fix a locus containing $Q$ variants indexed by $i in {1, dots, Q}$. For each trait $k in {1,2}$
there is a cohort of $N_k$ individuals with phenotype vector $y_k in RR^(N_k)$, additive dosage
vectors $g_(k i) in RR^(N_k)$, and a covariate matrix $X_k$ with $c_k$ columns, including an
intercept. The observed quantities are, for every $k$ and $i$, the coefficient $bh_(k i)$ and
standard error of the *marginal* regression of $y_k$ on $X_k$ and $g_(k i)$ alone. Write $V_(k i)$
for the squared standard error. We write $D_k$ for the data of cohort $k$ and $D = (D_1, D_2)$.
Everything is conditional on genotypes and covariates.

== Configurations and hypotheses

Let $S_k in {0, 1, dots, Q}$ be the *causal variant* of trait $k$ in the locus, with $S_k = 0$
meaning that the trait has none. Causal has a precise meaning here: $S_k = i$ says that $g_(k i)$
is the only variant in the locus with a nonzero coefficient in the joint regression of $y_k$. A
*configuration* is a pair $S = (S_1, S_2)$, and there are $(Q + 1)^2$ of them. They partition into
five sets, listed in @tab-hyp and drawn in @fig-grid.

#figure(
  table(
    columns: (auto, auto, auto, 1fr),
    align: (left, left, center, left),
    stroke: none,
    table.hline(),
    table.header([*Hypothesis*], [*Set of configurations*], [*Size $n_h$*], [*Meaning*]),
    table.hline(stroke: 0.5pt),
    [$H_0$], [$Sset_0 = {(0,0)}$], [$1$], [neither trait has a causal variant],
    [$H_1$], [$Sset_1 = {(i,0): i >= 1}$], [$Q$], [trait 1 only],
    [$H_2$], [$Sset_2 = {(0,j): j >= 1}$], [$Q$], [trait 2 only],
    [$H_3$], [$Sset_3 = {(i,j): i,j >= 1, thin i != j}$], [$Q(Q-1)$], [both traits, distinct variants],
    [$H_4$], [$Sset_4 = {(i,i): i >= 1}$], [$Q$], [both traits, one shared variant],
    table.hline(),
  ),
  caption: [The five hypotheses. The sizes sum to $(Q+1)^2$.],
) <tab-hyp>

#let cell-color(x, y) = {
  if x == 0 and y == 0 { rgb("#d9d9d9") } else if y == 0 { rgb("#f6d7a7") } else if x == 0 {
    rgb("#b7dcc8")
  } else if x == y { rgb("#d98b8b") } else { rgb("#bcd0e8") }
}
#figure(
  {
    let n = 6
    let labels = ([0], [1], [2], [3], [$dots$], [$Q$])
    let cells = ([$S_1 #h(2pt) \\ #h(2pt) S_2$],)
    for x in range(n) { cells.push(strong(labels.at(x))) }
    for y in range(n) {
      cells.push(strong(labels.at(y)))
      for x in range(n) {
        let name = if x == 0 and y == 0 { $H_0$ } else if x == 0 { $H_1$ } else if y == 0 {
          $H_2$
        } else if x == y { $H_4$ } else { $H_3$ }
        cells.push(table.cell(fill: cell-color(y, x), text(size: 8.5pt, name)))
      }
    }
    table(
      columns: (46pt,) + (30pt,) * n,
      rows: 19pt,
      align: center + horizon,
      stroke: (x, y) => if x > 0 and y > 0 { 0.4pt + white } else { none },
      ..cells,
    )
  },
  kind: image,
  caption: [The configuration space. Entry $(S_1, S_2)$ is colored by its hypothesis. The sum over
    $H_3$ is the sum over the whole inner square minus the diagonal, which is why it costs $O(Q)$
    rather than $O(Q^2)$ (@eq:h3sum).],
) <fig-grid>

== Likelihood

For a quantitative trait, given $S_k$ and an effect size $beta_k$,
$
  y_k = X_k alpha_k + g_(k S_k) beta_k + epsilon_k,
  quad epsilon_k dist Normal(0, sigma_k^2 I),
$ <eq:linear>
with the variant term absent when $S_k = 0$. For a case-control trait the same linear predictor
enters through a logit link, $op("logit") P(y_(k n) = 1) = x_(k n)^T alpha_k + g_(k n S_k) beta_k$.

== Assumptions

The whole method follows from five assumptions:

/ (A1) At most one causal variant: each trait has at most one causal variant in the locus. This is
  the encoding $S_k in {0, dots, Q}$.
/ (A2) Completeness: a causal variant, if there is one, is among the $Q$ variants analyzed.
/ (A3) Independent residuals: given $S$, the effects and the genotypes, $epsilon_1$ and $epsilon_2$
  are independent. Disjoint cohorts suffice.
/ (A4) Independent effects: given $S$, the effects $beta_1, beta_2$ are independent with
  $beta_k | {S_k != 0} dist Normal(0, omega_k^2)$. This holds even when $S_1 = S_2$.
/ (A5) Exchangeable variants: the prior on configurations is the one in @sec-prior, which treats
  all variants alike.

@sec-failures examines what happens when these fail. Giambartolomei et al. also assume that the
two cohorts share allele frequencies and LD [1]. The formulas below never use that: the dosage
vectors $g_(1 i)$ and $g_(2 i)$ are allowed to differ arbitrarily.

== Prior on configurations <sec-prior>

Following Text S1 of [1], give each variant, independently, one of four categories with
probabilities $p_0, p_1, p_2, p_(12)$ summing to 1: associated with neither trait, trait 1 only,
trait 2 only, or both. Let $Aev$ be the event that at most one variant is associated with each
trait. The prior on configurations is this product measure conditioned on $Aev$.

#lemma(title: "configuration odds")[
  Let $rho_1 = p_1 \/ p_0$, $rho_2 = p_2 \/ p_0$, $rho_(12) = p_(12) \/ p_0$, and let
  $S^0 = (0,0)$. For $S$ in $Sset_h$,
  $
    P(S | Aev) / P(S^0 | Aev) = cases(
      1 & h = 0,
      rho_1 & h = 1,
      rho_2 & h = 2,
      rho_1 rho_2 & h = 3,
      rho_(12) quad & h = 4.
    )
  $ <eq:configodds>
] <lem-prior>

#proof[
  A category vector satisfies $Aev$ exactly when at most one entry lies in {trait 1 only, both} and
  at most one in {trait 2 only, both}. Such vectors correspond one-to-one with configurations: $(i,j)$
  with $i != j$ is the vector with category "trait 1 only" at $i$ and "trait 2 only" at $j$, and
  $(i,i)$ is the vector with "both" at $i$. Their unconditional probabilities are $p_0^Q$,
  $p_0^(Q-1) p_1$, $p_0^(Q-1) p_2$, $p_0^(Q-2) p_1 p_2$ and $p_0^(Q-1) p_(12)$. Conditioning on $Aev$
  divides each by the same constant. Dividing by $p_0^Q$ gives @eq:configodds.
]

Consequently the prior on the hypotheses is, with $n_h$ from @tab-hyp and $rho_h$ the odds in
@eq:configodds,
$
  P(H_h | Aev) = (n_h rho_h) / Z,
  quad Z = 1 + Q rho_1 + Q rho_2 + Q(Q-1) rho_1 rho_2 + Q rho_(12).
$ <eq:hypprior>
The familiar approximations $P(H_1) approx Q p_1$ and $P(H_3) approx Q^2 p_1 p_2$ require
$Z approx 1$, which is not guaranteed. With the default $p_1 = p_2 = 10^(-4)$, $p_(12) = 10^(-5)$,
@eq:hypprior gives $(0.820, 0.082, 0.082, 0.008, 0.008)$ at $Q = 1000$ but
$(0.244, 0.244, 0.244, 0.244, 0.024)$ at $Q = 10^4$. The prior remains well defined for every $Q$,
because the conditional construction normalizes it. Only the reading of $Q p_1$ as a probability
breaks down.

= The Bayes factor for one trait <sec-abf>

Throughout this section there is one trait, so we drop $k$. The goal is the Bayes factor comparing
"variant $i$ is causal" against "no causal variant",
$
  BF(i) = p(D | S = i) / p(D | S = 0).
$

== A Gaussian integral

#lemma[
  For $a in RR$, $b >= 0$ and $omega > 0$,
  $
    integral_RR exp(a beta - (b beta^2) / 2) Normal(beta; 0, omega^2) dif beta
    = (1 + b omega^2)^(-1\/2) exp( (a^2 omega^2) / (2 (1 + b omega^2)) ).
  $ <eq:gauss>
] <lem-gauss>

#proof[
  Let $c = b + omega^(-2)$. The integrand is
  $(2 pi omega^2)^(-1\/2) exp(a beta - c beta^2 \/ 2)$. Completing the square,
  $a beta - c beta^2 \/ 2 = -c (beta - a\/c)^2 \/ 2 + a^2 \/ (2 c)$, so the integral equals
  $(2 pi omega^2)^(-1\/2) (2 pi \/ c)^(1\/2) exp(a^2 \/ (2c))$. Substituting
  $omega^2 c = 1 + b omega^2$ and $1\/c = omega^2 \/ (1 + b omega^2)$ gives @eq:gauss.
]

== Quantitative traits: the Bayes factor is exact <sec-exact>

Remove the covariates first. Let $U$ be a matrix whose $N - c$ rows are an orthonormal basis of the
orthogonal complement of the column space of $X$, and set $yt = U y$ and $gt_i = U g_i$. Applying
$U$ to @eq:linear eliminates $alpha$:
$
  yt | {S = i, beta} dist Normal(gt_i beta, sigma^2 I),
  quad yt | {S = 0} dist Normal(0, sigma^2 I).
$ <eq:projected>
By the Frisch-Waugh-Lovell theorem, the marginal regression of $y$ on $X$ and $g_i$ has coefficient
$
  bh_i = (gt_i^T yt) / norm(gt_i)^2,
  quad "with sampling variance" quad V_i = sigma^2 / norm(gt_i)^2.
$ <eq:ols>

#theorem(title: "the Bayes factor")[
  Under @eq:projected with $sigma^2$ known and $beta dist Normal(0, omega^2)$, for every genotype
  matrix, regardless of LD,
  $
    BF(i) = sqrt(1 - r_i) thin exp((r_i z_i^2) / 2),
    quad r_i = omega^2 / (omega^2 + V_i),
    quad z_i = bh_i / sqrt(V_i).
  $ <eq:abf>
] <thm-abf>

#proof[
  For fixed $beta$, the ratio of the densities in @eq:projected is
  $
    exp( - (norm(yt - gt_i beta)^2 - norm(yt)^2) / (2 sigma^2) )
    = exp( (gt_i^T yt) / sigma^2 thin beta - norm(gt_i)^2 / (2 sigma^2) thin beta^2 )
    = exp( bh_i / V_i thin beta - beta^2 / (2 V_i) ),
  $ <eq:lrform>
  using @eq:ols in the last step. $BF(i)$ is the average of @eq:lrform over the prior, which is the
  integral of @lem-gauss with $a = bh_i \/ V_i$ and $b = 1 \/ V_i$. Now
  $1 + b omega^2 = (V_i + omega^2) \/ V_i = 1 \/ (1 - r_i)$, and since
  $omega^2 (1 - r_i) = r_i V_i$,
  $
    (a^2 omega^2) / (2 (1 + b omega^2))
    = bh_i^2 / (2 V_i^2) thin omega^2 (1 - r_i)
    = bh_i^2 / (2 V_i^2) thin r_i V_i
    = (r_i z_i^2) / 2.
  $
]

We call @eq:abf the approximate Bayes factor $ABF(i)$, following the literature, and write
$Lambda_i = log ABF(i) = 1/2 [log(1 - r_i) + r_i z_i^2]$. The name is slightly misleading: under
the model of this section there is nothing approximate about it. The approximations are the two
that follow.

#remark[why LD disappears][
  This is the point of the method, and it deserves emphasis. The exponent in @eq:lrform involves
  $y$ and $g_i$ but no other column. Hypothesis $S = i$ differs from the common reference $S = 0$
  only in the direction $gt_i$, so the likelihood ratio sees only the projection of the data on that
  direction. LD is not ignored by an approximation. It is genuinely absent from each ratio against
  the null. LD still shapes the *values* of the $z_i$ (@sec-ld), but it is not needed to convert
  those values into Bayes factors.

  This depends entirely on (A1). If two variants $i, j$ are jointly causal, the same computation
  gives
  $
    log (p(yt | beta_i, beta_j)) / (p(yt | 0))
    = (beta_i gt_i^T yt + beta_j gt_j^T yt) / sigma^2
    - (beta_i^2 norm(gt_i)^2 + beta_j^2 norm(gt_j)^2 + 2 beta_i beta_j thin gt_i^T gt_j)
    / (2 sigma^2),
  $
  and the cross term $gt_i^T gt_j$ is LD. Allowing multiple causal variants therefore requires an LD
  matrix. That is the real price of dropping (A1).
]

#remark[sufficiency][
  An equivalent proof explains why a single number suffices. Split $yt$ into its component along
  $gt_i$, which is $bh_i gt_i$, and the orthogonal remainder. The two are independent, and the
  remainder has the same law under $S = i$ and $S = 0$. So it cancels in the ratio, leaving
  $
    BF(i) = (Normal(bh_i; 0, V_i + omega^2)) / (Normal(bh_i; 0, V_i)).
  $ <eq:abfratio>
  This uses $bh_i | S = i dist Normal(0, V_i + omega^2)$ marginally over $beta$, and expanding the
  ratio recovers @eq:abf. Note what is *not* claimed. The $bh_i$ are correlated across variants, but
  each $BF(i)$ needs only the one $bh_i$, so that correlation never enters. The posterior of the
  effect under $S = i$ is $beta | yt dist Normal(r_i bh_i, r_i V_i)$, so $r_i$ is a shrinkage factor.
]

#remark[covariate prior][
  The reduction to @eq:projected is the same as placing a flat prior on $alpha$. That integral
  factors into a term in $yt$ and a constant that does not depend on $beta$, $i$ or the model, so
  the improper prior causes no ambiguity in the Bayes factor.
]

#remark[first approximation: plug-in variance][
  In practice $sigma^2$ is unknown and $V_i$ is the reported squared standard error, computed with
  a residual variance estimated in the model containing $g_i$. Let $e_i$ be the fraction of residual
  phenotypic variance explained by variant $i$, so that
  $z_i^2 = (N - c - 1) e_i \/ (1 - e_i)$. If $sigma^2$ were profiled out, twice the log likelihood
  ratio would be $-(N - c) log(1 - e_i)$. The two agree to first order in $e_i$, but
  $e_i\/(1-e_i)$ exceeds $-log(1-e_i)$: they are $0.111$ versus $0.105$ at $e_i = 0.1$ and $0.429$
  versus $0.357$ at $e_i = 0.3$. So the plug-in slightly overstates the strongest single-variant
  signals. It is negligible except for variants explaining a substantial share of the trait, and it
  preserves the ordering of variants by $e_i$.
]

== Case-control traits: the Bayes factor is asymptotic <sec-logistic>

For the logistic model there is no exact reduction to $(bh_i, V_i)$. The structural point survives:
under $S = i$ the likelihood involves $g_i$ and no other variant, and the reference $S = 0$ is
common to all $i$. The Gaussian form returns through a local quadratic expansion. Let
$ell(alpha, beta)$ be the log likelihood of the regression on $X$ and $g_i$, let
$tilde(ell)(beta) = max_alpha ell(alpha, beta)$ be its profile, and let $bh$ be the maximizer. The
Wald variance is $V = -1 \/ tilde(ell)''(bh)$. This is the $(beta, beta)$ entry of the inverse
information: by the Schur complement, profiling out $alpha$ yields the same curvature as the
reparameterization that makes the intercept and $beta$ orthogonal, used in Text S1 of [1]. Taylor
expansion gives
$
  tilde(ell)(beta) = tilde(ell)(bh) - (beta - bh)^2 / (2 V) + R(beta),
  quad R(beta) = O(N abs(beta - bh)^3).
$
Subtracting the value at $beta = 0$,
$
  tilde(ell)(beta) - tilde(ell)(0) = bh / V thin beta - beta^2 / (2 V) + R(beta) - R(0).
$ <eq:logisticquad>
Up to the remainder, this is exactly the exponent of @eq:lrform. A Laplace approximation of the
integral over $alpha$ with a flat prior contributes a further factor $1 + O(abs(beta))$ from the
change in curvature. Averaging over $beta dist Normal(0, omega^2)$ with @lem-gauss then gives
@eq:abf again.

The error is controlled as follows. The posterior of $beta$ concentrates within $O(sqrt(V))$ of
$r bh$, where $R(beta) = O(N V^(3\/2)) = O(N^(-1\/2))$, while $R(0) = O(N abs(bh)^3) =
O(z^2 abs(bh))$. So
$
  log BF(i) = Lambda_i + O(z_i^2 abs(bh_i)) + O(abs(bh_i)) + O(N^(-1\/2)).
$
For a fixed $z_i$ the error vanishes as $N -> oo$, and for fixed $bh_i$ it is a small relative error
in $z_i^2$. A direct check by nested numerical integration of the exact logistic likelihood, with
$N = 20000$, allele frequency $0.3$, about 30% cases and $omega = 0.2$:

#figure(
  table(
    columns: 5,
    align: right,
    stroke: none,
    table.hline(),
    table.header([*true $beta$*], [*$z$*], [*$Lambda$ (ABF)*], [*exact $log BF$*], [*difference*]),
    table.hline(stroke: 0.5pt),
    [0.00], [1.98], [$-0.215$], [$-0.219$], [0.004],
    [0.05], [$-0.38$], [$-2.066$], [$-2.066$], [0.000],
    [0.10], [2.64], [1.304], [1.295], [0.009],
    [0.20], [9.23], [39.904], [39.658], [0.246],
    [0.40], [16.50], [132.072], [131.803], [0.269],
    table.hline(),
  ),
  caption: [Accuracy of the approximate Bayes factor for logistic regression (one simulated data set per row).],
) <tab-logistic>

Case-control ascertainment does not disturb this. Under retrospective sampling the logistic slope
keeps its meaning and its estimator keeps its asymptotic distribution, and only the intercept
shifts (Prentice and Pyke [6]).

== Reading the Bayes factor <sec-reading>

Since $(1 - r_i) z_i^2 = bh_i^2 \/ (V_i + omega^2)$, the log Bayes factor separates into three terms:
$
  Lambda_i = underbrace(1/2 log V_i / (V_i + omega^2), "Occam factor") +
  underbrace(z_i^2 / 2, "likelihood gain") -
  underbrace(bh_i^2 / (2 (V_i + omega^2)), "prior-tail penalty").
$ <eq:decomp>
The first term is negative and grows in magnitude as precision improves, so a variant with
$z_i = 0$ counts against association: the data support the null, and more so in larger samples. The
second is the usual evidence. The third penalizes estimates far in the tail of the prior. It is
appreciable only when $abs(bh_i)$ is comparable to $omega$ or larger.

The Bayes factor exceeds 1 exactly when $z_i^2 > -log(1 - r_i) \/ r_i$. For $sqrt(V_i) = 0.05$ and
$omega = 0.15$ that is $abs(z_i) > 1.60$, and for $sqrt(V_i) = 0.01$ it is $abs(z_i) > 2.33$.
Because of $r_i$ and the first term, $Lambda_i$ is not a function of $z_i$ alone: among variants
with equal $abs(z_i)$, the Bayes factors differ according to precision.

= The posterior over hypotheses <sec-posterior>

Write $B_k (i) = ABF_k (i)$ for $i >= 1$, with the convention $B_k (0) = 1$, and define
$
  L_1 = sum_(i=1)^Q B_1 (i),
  quad L_2 = sum_(i=1)^Q B_2 (i),
  quad L_(12) = sum_(i=1)^Q B_1 (i) B_2 (i).
$ <eq:ls>

#theorem(title: "posterior probabilities")[
  Under (A1) to (A5), with $B_k$ the Bayes factors of @sec-abf,
  $
    P(H_h | D, Aev) = w_h / (sum_(ell=0)^4 w_ell),
  $ <eq:posterior>
  $
    w_0 = 1, quad w_1 = rho_1 L_1, quad w_2 = rho_2 L_2,
    quad w_3 = rho_1 rho_2 (L_1 L_2 - L_(12)),
    quad w_4 = rho_(12) L_(12).
  $ <eq:weights>
] <thm-posterior>

#proof[
  By Bayes' theorem, dividing the numerator and denominator by $P(S^0 | Aev) thin p(D | S^0)$,
  $
    P(S | D, Aev) prop P(S | Aev) / P(S^0 | Aev) thin dot.c thin p(D | S) / p(D | S^0).
  $ <eq:posteriorconfig>
  This is Eq. 3 of [1]. The first factor is @eq:configodds. For the second, integrate over the
  effects. By (A3) the phenotypes are conditionally independent given the effects, and by (A4) the
  effects have a product prior, so the double integral factors, including when $S_1 = S_2$:
  $
    p(D | S)
    &= integral.double p(D_1 | S_1, beta_1) p(D_2 | S_2, beta_2)
    thin pi(beta_1 | S_1) pi(beta_2 | S_2) dif beta_1 dif beta_2 \
    &= p(D_1 | S_1) thin p(D_2 | S_2).
  $
  Hence $p(D | S) \/ p(D | S^0) = B_1 (S_1) B_2 (S_2)$, using $B_k (0) = 1$. Summing
  @eq:posteriorconfig over each $Sset_h$ gives weight $1$ on $Sset_0$, then $rho_1 sum_i B_1 (i)$
  and $rho_2 sum_j B_2 (j)$ on $Sset_1$ and $Sset_2$, and $rho_(12) sum_i B_1 (i) B_2 (i)$ on
  $Sset_4$. On $Sset_3$ the sum runs over ordered pairs of distinct variants:
  $
    sum_(i != j) B_1 (i) B_2 (j)
    = sum_(i, j) B_1 (i) B_2 (j) - sum_i B_1 (i) B_2 (i)
    = L_1 L_2 - L_(12).
  $ <eq:h3sum>
  Normalizing over the five hypotheses gives @eq:posterior.
]

The complete $(Q+1)^2$ enumeration thus reduces to three $O(Q)$ sums. The weight $w_3$ is
nonnegative because it is a sum of positive terms, and it vanishes only when $Q = 1$.

#corollary(title: "variant-level posteriors")[
  Given $H_4$, the shared variant has posterior
  $P(S_1 = S_2 = i | H_4, D) = B_1 (i) B_2 (i) \/ L_(12)$. Given $H_3$,
  $P(S_1 = i, S_2 = j | H_3, D) = B_1 (i) B_2 (j) \/ (L_1 L_2 - L_(12))$ for $i != j$.
]

#remark[sign blindness][
  $Lambda_(k i)$ depends on $bh_(k i)$ only through its square. So the posterior is invariant to
  flipping the sign of any single $bh_(k i)$. The method does not use whether the two traits' effects
  have consistent directions across a locus, and allele orientation matters only for deciding which
  rows are the same variant. This is a direct consequence of the centered product prior (A4).
]

== The prior actually used in implementations <sec-implemented>

The R package and `coloc_abf.py` both use the weights
$(1, p_1 L_1, p_2 L_2, p_1 p_2 (L_1 L_2 - L_(12)), p_(12) L_(12))$, that is, @eq:weights with $p$ in
place of $rho$. Text S1 of [1] justifies this as the approximation $p_0 approx 1$. A sharper
statement is available.

#proposition(title: "the implemented posterior is exact for a rescaled prior")[
  The weights above are exactly @eq:weights for the per-variant probabilities
  $p'_bullet = p_bullet \/ (1 + p_1 + p_2 + p_(12))$ and
  $p'_0 = 1 \/ (1 + p_1 + p_2 + p_(12))$.
] <prop-implemented>

#proof[
  These four numbers sum to 1, and $p'_bullet \/ p'_0 = p_bullet$. So $rho_1 = p_1$, $rho_2 = p_2$,
  $rho_(12) = p_(12)$, and hence $rho_1 rho_2 = p_1 p_2$.
]

So the implemented posterior is an exact posterior, and $p_1, p_2, p_(12)$ are best read as
*per-variant odds against the "neither" category*. At the defaults $p_0 = 0.99979$, so the
distinction is numerically immaterial. It is worth stating because it shows that no constraint on
$Q$ is needed.

#remark[large regions][
  The R package replaces any $p$ with $1 \/ (Q + 1)$ when $Q p >= 1$, a guard it calls
  `adjust_prior`. `coloc_abf.py` has no such guard. At the defaults the two first differ at
  $Q = 10^4$ (for $p_1, p_2$). @prop-implemented shows the unadjusted computation is still a valid
  posterior there, so this is a difference of convention rather than an error in either.
]

#remark[values of the priors][
  The main analysis of [1] uses $p_(12) = 10^(-6)$, with $2 times 10^(-6)$ and $10^(-5)$ as
  sensitivity analyses. The R package and `coloc_abf.py` default to $10^(-5)$. Wallace [3], using
  simulation, reports that $10^(-5)$ may be overly liberal and suggests $5 times 10^(-6)$.
  @sec-odds quantifies what these choices mean.
]

= Inputs: sampling variance and prior scale <sec-inputs>

@eq:abf needs $V_i$ and $omega$. The reported standard error supplies $V_i$ directly. What
remains is to choose $omega$, and to approximate $V_i$ where it is not reported.

We use $omega$ for the prior standard deviation. This is $W$ in `coloc_abf.py`, whereas in [1] the
symbol $W$ denotes the prior *variance* $omega^2$. Text S1 of [1] calls the value $0.2$ a variance
in one sentence, but its main text and the R code both treat $0.15$ and $0.2$ as standard
deviations.

== Quantitative traits

The prior scale must carry the units of the phenotype.

#proposition(title: "scale invariance")[
  Under $y -> c y$, the Bayes factor @eq:abf is unchanged for all data if and only if
  $omega -> abs(c) omega$.
] <prop-scale>

#proof[
  $bh_i -> c bh_i$ and $V_i -> c^2 V_i$, so $z_i^2$ is invariant. Thus $Lambda_i$ changes only through
  $r_i$. For $Lambda_i = 1/2[log(1-r_i) + r_i z_i^2]$ to be unchanged for all values of $z_i^2$, the
  coefficient $r_i$ of $z_i^2$ must be unchanged, which means $omega^2 \/ V_i$ is unchanged, which
  means $omega -> abs(c) omega$. The converse is immediate.
]

So set $omega = 0.15 thin "sd"(y)$. Under this prior, 95% of causal effects lie within $plus.minus
0.29$ phenotypic standard deviations. An effect of exactly $omega$ at allele frequency $f$ explains a
fraction $omega^2 dot.c 2 f (1-f) \/ "sd"(y)^2$ of phenotypic variance, about $0.01$ at $f = 0.3$.

When $"sd"(y)$ is unknown it can be recovered from the summary statistics. From @eq:ols, if the
variant explains little variance then $sigma^2 approx "sd"(y)^2$, and under Hardy-Weinberg
equilibrium $norm(gt_i)^2 approx N_i dot.c 2 f_i (1 - f_i)$, because $g_i dist "Binomial"(2, f_i)$.
So
$
  V_i approx "sd"(y)^2 / (2 N_i f_i (1 - f_i)),
  quad "that is," quad
  underbrace(2 N_i f_i (1 - f_i), v_i) approx "sd"(y)^2 dot.c underbrace(1 / V_i, u_i).
$ <eq:vquant>
This is a line through the origin with slope $"sd"(y)^2$. Least squares minimizes
$sum_i (v_i - c thin u_i)^2$, which gives
$
  hat("sd")(y)^2 = (sum_i u_i v_i) / (sum_i u_i^2).
$ <eq:sdy>
The estimate is positive, because all $u_i, v_i$ are. The formula uses $f_i (1-f_i)$, which is
symmetric, so an effect allele frequency works as well as a minor allele frequency.

This estimate is not delicate. Differentiating $Lambda_i$,
$
  (partial Lambda_i) / (partial log omega) = -r_i + r_i (1 - r_i) z_i^2,
  quad "which is about" -1 + bh_i^2 / omega^2 "when" V_i << omega^2.
$
So, among variants with $abs(bh_i)$ well below $omega$ or with similar $abs(bh_i)$, a small relative
error in $omega$ shifts all $Lambda_i$ by nearly the same amount. That common shift changes $L_k$ but
leaves the variant-level posteriors $pi_k$ of @sec-reduction essentially unchanged.

== Case-control traits

Here $beta$ is a log odds ratio, which has no units, and [1] takes $omega = 0.2$. Under this prior
95% of odds ratios lie between $0.68$ and $1.48$.

If $V_i$ is not reported, it can be approximated from the information matrix at the null. With
success probability $s$ (the case fraction) for every individual, the entries for the intercept and
$beta$ are $I_(alpha alpha) = N s(1-s)$, $I_(alpha beta) = s(1-s) sum_n g_n$ and
$I_(beta beta) = s(1-s) sum_n g_n^2$. The Schur complement gives
$
  1 / V_i = I_(beta beta) - I_(alpha beta)^2 / I_(alpha alpha)
  = s(1-s) sum_n (g_n - overline(g))^2
  approx N s (1 - s) dot.c 2 f_i (1 - f_i).
$ <eq:vcc>
So $s$ enters only through $V_i$. When standard errors are reported, $s$ plays no part in the
calculation. Accordingly, `coloc_abf.py` stores `case_fraction` without using it.

= The structure of the posterior <sec-structure>

== Reduction to four numbers <sec-reduction>

Define, for $k = 1, 2$,
$
  pi_k (i) = (B_k (i)) / L_k,
  quad a_k = (rho_k L_k) / (1 + rho_k L_k),
  quad C = sum_(i=1)^Q pi_1 (i) pi_2 (i) = L_(12) / (L_1 L_2),
  quad kappa = rho_(12) / (rho_1 rho_2).
$ <eq:reduction>
Each has a plain meaning:

- $pi_k$ is the *fine-mapping posterior* of trait $k$: the posterior on the location of its causal
  variant, given that it has one.
- $a_k$ is the probability that trait $k$ has a causal variant in the locus, according to that
  trait's data alone, with per-variant prior odds $rho_k$.
- $C$ is the *overlap*: the probability that independent draws from $pi_1$ and $pi_2$ coincide.
- $kappa = p_(12) p_0 \/ (p_1 p_2)$ is the *enrichment*: the odds ratio of the per-variant $2 times
  2$ table of (associated with trait 1) against (associated with trait 2). $kappa = 1$ means that
  the two kinds of association fall on variants independently, and the default is $kappa = 1000$.

#theorem(title: "four-number reduction")[
  Let $T = 1 + (kappa - 1) a_1 a_2 C$. Then
  $
    PP_0 = ((1 - a_1)(1 - a_2)) / T, quad
    PP_1 = (a_1 (1 - a_2)) / T, quad
    PP_2 = ((1 - a_1) a_2) / T,
  $
  $
    PP_3 = (a_1 a_2 (1 - C)) / T, quad
    PP_4 = (kappa a_1 a_2 C) / T.
  $ <eq:reduced>
] <thm-reduction>

#proof[
  Expanding the sum of the weights in @eq:weights,
  $
    sum_(h=0)^4 w_h
    &= 1 + rho_1 L_1 + rho_2 L_2 + rho_1 rho_2 L_1 L_2 + (rho_(12) - rho_1 rho_2) L_(12) \
    &= (1 + rho_1 L_1)(1 + rho_2 L_2) + (rho_(12) - rho_1 rho_2) L_(12).
  $ <eq:factor>
  Divide every $w_h$ and @eq:factor by $(1 + rho_1 L_1)(1 + rho_2 L_2)$. Since
  $1 \/ (1 + rho_k L_k) = 1 - a_k$ and
  $rho_1 rho_2 L_1 L_2 \/ ((1 + rho_1 L_1)(1 + rho_2 L_2)) = a_1 a_2$, while
  $L_(12) = C L_1 L_2$ and $rho_(12) = kappa rho_1 rho_2$, the denominators become $T$ and the
  numerators become those of @eq:reduced.
]

So the posterior depends on the data through only $a_1$, $a_2$ and $C$. When $kappa = 1$, $T = 1$
and the posterior is a product of two independent single-trait analyses, with $PP_4 = a_1 a_2 C$
being the chance that two unrelated signals happen to land on one variant. The role of $kappa$ is
to multiply that chance.

== The odds of $H_4$ against $H_3$ <sec-odds>

When both traits clearly have a signal, the question is $H_4$ against $H_3$, and @eq:reduced gives
a very simple answer.

#corollary(title: "the central ratio")[
  $
    PP_4 / PP_3 = kappa thin C / (1 - C).
  $ <eq:ratio>
  The prior odds are $kappa \/ (Q - 1)$, so the Bayes factor of $H_4$ against $H_3$ is
  $(Q - 1) C \/ (1 - C)$. It equals 1 at the uniform baseline $C = 1 \/ Q$.
] <cor-ratio>

Thus the entire comparison of a shared variant against distinct variants depends on the data only
through the overlap of the two fine-mapping posteriors. Three consequences follow.

*How much overlap is enough.* Since $PP_4 \/ (PP_3 + PP_4) >= tau$ exactly when
$
  C >= tau / (tau + kappa (1 - tau)),
$ <eq:cthreshold>
the bar is low. For $tau = 0.8$ the threshold is $C >= 0.0040$ at $p_(12) = 10^(-5)$, $0.0079$ at
$5 times 10^(-6)$, and $0.038$ at $10^(-6)$. With $Q = 2000$, the default therefore needs only
eight times the overlap expected from two uniform posteriors. This makes precise the concern in [3]
that the default may be liberal: $H_4$ is favored unless the data show that the two signals are
*separated*.

*An upper bound from resolution.* By Cauchy-Schwarz, $C <= (n_1 n_2)^(-1\/2)$, where
$n_k = 1 \/ sum_i pi_k (i)^2$ is the effective number of candidate variants for trait $k$. So
$PP_4 \/ PP_3 <= kappa \/ (sqrt(n_1 n_2) - 1)$. Poorly resolved signals cap the attainable support
for $H_4$.

*The number of variants has no explicit role.* $Q$ appears in neither the weights nor @eq:ratio. The
prior odds fall like $1 \/ (Q - 1)$ (Wallace [3]), but the Bayes factor rises by exactly the
compensating factor. Adding variants with no signal changes the $L$'s only slightly, so it leaves
$pi_k$ nearly unchanged and the posterior odds of $H_4$ against $H_3$ with them. The parameter that
really matters is the per-variant $kappa$, not the region-level prior.

== Sensitivity to $p_(12)$ and borrowing of strength

#corollary(title: "sensitivity")[
  $p_(12)$ enters only $w_4$, so $PP_4 \/ (1 - PP_4)$ is proportional to $rho_(12)$. From a single run
  at $rho_(12)$, the value of $PP_4$ at any other $rho'_(12)$ follows from
  $
    (PP_4 (rho'_(12))) / (1 - PP_4 (rho'_(12)))
    = rho'_(12) / rho_(12) dot.c (PP_4 (rho_(12))) / (1 - PP_4 (rho_(12))).
  $
  This is the single-parameter case of the sensitivity reweighting in [3].
]

Large $kappa$ also lets a strong signal support a weak one. If $a_1 approx 1$, then
$PP_4 \/ PP_1 approx kappa a_2 C \/ (1 - a_2)$, so $a_2$ of order $1 \/ (kappa C)$ is enough for $H_4$
to overtake $H_1$. A trait with no significant signal on its own can therefore be assigned to a
shared variant. The "trait 2 weak" row of @tab-stylized shows this: with $a_2 = 0.035$, the
posterior is still $PP_4 = 0.73$.

== How LD acts through the $z$ statistics <sec-ld>

LD is not an input, but it determines the distribution of what is observed. Let $R_(i j)$ be the
correlation of $gt_i$ and $gt_j$, the sample LD.

#proposition(title: [law of the $z$ statistics])[
  Under @eq:projected with $S = c$, and with $z_i = gt_i^T yt \/ (sigma norm(gt_i))$,
  $
    (z_1, dots, z_Q) dist Normal(lambda R_(dot.c c), R),
    quad lambda = beta / sqrt(V_c).
  $ <eq:zlaw>
] <prop-zlaw>

#proof[
  With unit vectors $t_i = gt_i \/ norm(gt_i)$ we have $z_i = t_i^T yt \/ sigma$ and
  $yt dist Normal(beta gt_c, sigma^2 I)$. So
  $EE z_i = beta norm(gt_c) t_i^T t_c \/ sigma = lambda R_(i c)$ and
  $"Cov"(z_i, z_j) = t_i^T t_j = R_(i j)$.
]

Therefore $EE[z_c^2 - z_i^2] = lambda^2 (1 - R_(i c)^2)$. When $r_i approx 1$ for all variants,
$log(pi_k (c) \/ pi_k (i)) approx (z_c^2 - z_i^2) \/ 2$, so a variant in LD $R_(i c)^2$ with the
causal variant remains a live candidate as long as $lambda^2 (1 - R_(i c)^2)$ is of order 1. This
is how LD spreads $pi_k$, and it controls the attainable $C$.

#proposition(title: "perfect LD")[
  Suppose a set $M$ of $m$ variants has $R_(i j)^2 = 1$ and equal $V_i$ within each trait, and that
  $pi_1$ and $pi_2$ are concentrated on $M$. Then $C = 1 \/ m$ and
  $PP_4 \/ PP_3 = kappa \/ (m - 1)$. The Bayes factor of $H_4$ against $H_3$ is $(Q-1)\/(m-1)$.
] <prop-clique>

#proof[
  Perfect correlation gives $z_i^2$ equal on $M$, and with equal $V_i$ the Bayes factors are equal
  there. So $pi_k$ is uniform on $M$, and $C = m dot.c m^(-2)$. Substitute in @cor-ratio.
]

Here the data cannot distinguish a shared variant from two distinct variants inside $M$, and the
result is a prior calculation restricted to $M$. With $m - 1 < kappa$ the default prior chooses
$H_4$. So $PP_4$ is better read as "the two signals are not separable" than as a proof of a single
shared molecular variant.

#proposition(title: "consistency")[
  Fix $Q$ and the true configuration, with fixed nonzero effects. Let $N_1, N_2 -> oo$ with
  $V_(k i) = v_(k i) \/ N_k$ for fixed $v_(k i)$ and the sample LD converging to a limit $R$. If
  $R_(i c)^2 < 1$ for each causal variant $c$ and every $i != c$, then the posterior probability of
  the true hypothesis tends to 1 in probability.
] <prop-consistency>

#proof[
  If trait $k$ has no causal variant, then $z_(k i) = O_p (1)$ by @prop-zlaw and
  $1 - r_(k i) tilde.op V_(k i) \/ omega_k^2$, so $B_k (i) = O_p (N_k^(-1\/2))$. Thus $L_k -> 0$ and
  $a_k -> 0$. If trait $k$ has causal variant $c$, then by @prop-zlaw
  $z_(k i)^2 \/ N_k -> delta_k R_(i c)^2$ with $delta_k = beta_k^2 \/ v_(k c) > 0$. So
  $Lambda_(k c) - Lambda_(k i) -> oo$ for $i != c$ and $Lambda_(k c) -> oo$, each linearly in $N_k$
  up to logarithmic terms, which gives
  $pi_k -> bold(e)_c$ and $a_k -> 1$. Substitute in @eq:reduced. For truth $H_0$, $H_1$ or $H_2$, one
  of the $a_k$ tends to zero, so $T -> 1$ and the claim is immediate. For truth $H_3$, $C -> 0$, so
  $T -> 1$ and $PP_3 -> 1$. For truth $H_4$, $C -> 1$, so $T -> kappa$ and $PP_4 -> 1$.
]

If the condition fails because of a clique of size $m$ as in @prop-clique, the same argument gives
the limit $PP_4 -> kappa \/ (kappa + m - 1)$ under true $H_4$. So the posterior does not concentrate,
and no amount of data helps.

== Stylized examples

@tab-stylized applies @eq:reduced to noise-free data: $z_(k i)$ is set to its expectation from
@prop-zlaw, with $Q = 1000$, $R_(i j) = 0.98^(abs(i - j))$, $sqrt(V) = 0.01$, $omega = 0.15$,
default priors and peak $abs(z) = 10$ unless stated.

#figure(
  table(
    columns: (1fr, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right),
    stroke: none,
    table.hline(),
    table.header([*Locus*], [$a_1$], [$a_2$], [$C$], [$PP_1$], [$PP_3$], [$PP_4$]),
    table.hline(stroke: 0.5pt),
    [Shared variant], [1.000], [1.000], [0.589], [0.000], [0.001], [0.999],
    [Distinct variants, 100 positions apart], [1.000], [1.000], [$8.7 times 10^(-22)$], [0.000], [1.000], [0.000],
    [Distinct variants, adjacent ($R = 0.98$)], [1.000], [1.000], [0.161], [0.000], [0.005], [0.995],
    [Trait 2 null], [1.000], [0.007], [0.001], [0.987], [0.007], [0.007],
    [Shared variant, trait 2 weak ($z = 3.5$)], [1.000], [0.035], [0.076], [0.262], [0.009], [0.729],
    [Trait 1 has two signals ($z = 10, 8$); trait 2 shares the weaker], [1.000], [1.000], [$1.0 times 10^(-8)$], [0.000], [1.000], [0.000],
    table.hline(),
  ),
  caption: [Posterior on noise-free loci. $PP_0$ and $PP_2$ are 0.000 throughout. The third row has
    truth $H_3$ but $PP_4 = 0.995$, and the sixth has a truly shared signal but $PP_3 = 1$.],
) <tab-stylized>

The third and sixth rows are the two characteristic errors. The third follows from @eq:ratio:
$C = 0.16$ is far above the threshold @eq:cthreshold, and these signals cannot be separated at this
power. The sixth is the subject of the next section.

= When the assumptions fail <sec-failures>

== More than one causal variant (A1)

#proposition(title: "the masked signal")[
  Suppose trait 1 has two signals at unlinked variants $a, b$ with $abs(z_(1 a)) > abs(z_(1 b))$, and
  that $pi_1$ is concentrated on ${a, b}$ with $r$ common to both. Suppose that $pi_2$ is concentrated
  on $b$ alone. Then
  $
    PP_4 / PP_3 = kappa exp( - r / 2 (z_(1 a)^2 - z_(1 b)^2) ).
  $ <eq:masked>
  So $H_3$ is preferred once $z_(1 a)^2 - z_(1 b)^2 > 2 log(kappa) \/ r$, which at the default is
  about $13.8$.
] <prop-masked>

#proof[
  Here $pi_1 (b) \/ pi_1 (a) = exp(-r(z_(1 a)^2 - z_(1 b)^2) \/ 2)$ and $C = pi_1 (b)$, so
  $C \/ (1 - C) = pi_1 (b) \/ pi_1 (a)$. Apply @cor-ratio.
]

For $(z_(1 a), z_(1 b)) = (9, 8)$ this gives $PP_4 \/ PP_3 = 0.21$, and for $(10, 8)$ it gives
$1.6 times 10^(-5)$. The weaker signal is masked, because (A1) forces the posterior to place the
one causal variant of trait 1 on the stronger peak. coloc then reports distinct variants even
though $b$ really is shared. The report is correct about the stronger signal of trait 1, and what is
lost is the shared weaker signal. This is intrinsic to the method, and by the first remark of
@sec-exact it cannot be fixed without LD.

== Correlated residuals (A3)

For individuals measured on both traits with residual correlation $rho_e$, the joint log likelihood
ratio for configuration $(i, j)$ contains the term
$
  (rho_e) / (sigma_1 sigma_2 (1 - rho_e^2)) thin beta_1 beta_2 sum_n g_(n i) g_(n j).
$
So the factorization in the proof of @thm-posterior fails, and LD between $i$ and $j$ re-enters. The
precise condition is on the residuals. Sharing individuals does no harm if $rho_e = 0$, because
everything is conditional on genotypes. Conversely, disjoint cohorts are a sufficient condition,
which is the setting of [1]. Shared controls in two case-control studies would induce $rho_e != 0$.

== A missing causal variant (A2)

If the true causal variant $c$ is not among the $Q$ analyzed, then @eq:linear does not hold for any
$S$. By the argument of @prop-zlaw, each trait's $z$ statistics still have mean proportional to the
correlation with $c$, so with equal LD in the two cohorts the two peaks tend to fall on the same
tagging variants, and $H_4$ can still be supported. Sparse or unequal variant coverage is the main
risk. The simulations of [1] (Figures S1 and S2) examine this case.

== Independent effects (A4)

A real shared causal variant may have correlated effects on the two traits. The product prior does
not use this. That is conservative rather than invalid: it forgoes information from concordance of
effect sizes (see the sign-blindness remark in @sec-posterior).

= Numerical evaluation <sec-numerics>

Individual Bayes factors can exceed the range of floating point numbers, since $z^2 \/ 2$ can reach
the thousands. So every step is carried out on the log scale. Let
$"lse"(x) = m + log sum_i e^(x_i - m)$ with $m = max_i x_i$. Then:

+ For each trait, compute $r_(k i)$, $z_(k i)$ and $Lambda_(k i)$ from @eq:abf. For numerical
  accuracy, use $log(1 - r) = "log1p"(-r)$.
+ Compute $ell_1 = "lse"(Lambda_1)$, $ell_2 = "lse"(Lambda_2)$ and $ell_(12) = "lse"(Lambda_1 + Lambda_2)$.
+ Set $log w_0 = 0$, $log w_1 = log p_1 + ell_1$, $log w_2 = log p_2 + ell_2$ and
  $log w_4 = log p_(12) + ell_(12)$.
+ With $a = ell_1 + ell_2$ and $b = ell_(12)$, set
  $log w_3 = log p_1 + log p_2 + a + "log1p"(-e^(b - a))$, or $-oo$ if $b >= a$.
+ Return $exp(log w_h - "lse"(log w))$.

@eq:h3sum guarantees $b <= a$, with equality exactly when $Q = 1$. So the case $b >= a$ arises
only for $Q = 1$ or from rounding when $1 - C$ is below machine resolution relative to $abs(a)$.
Setting $w_3 = 0$ then is harmless, since by @eq:ratio $PP_3 \/ PP_4 = (1 - C) \/ (kappa C)$ is itself
negligible. The closed forms of @thm-reduction are for analysis. Computing $1 - C$ directly from
exponentiated logs would lose this precision.

#pagebreak(weak: true)
#heading(numbering: none)[Sources and verification]

*Sources read for this note.*

+ Giambartolomei C, Vukcevic D, Schadt EE, Franke L, Hingorani AD, Wallace C, Plagnol V. Bayesian test
  for colocalisation between pairs of genetic association studies using summary statistics.
  _PLoS Genetics_ 10(5):e1004383, 2014. Main text and Text S1. Source of the hypotheses, the per-variant prior,
  Eq. 3 (our @eq:posteriorconfig), and the values $omega = 0.15$ and $0.2$.
+ The `coloc` R package, file `R/claudia.R` on the main branch of `chr1swallace/coloc`: functions
  `approx.bf.estimates`, `combine.abf`, `sdY.est`, `Var.data.cc`, `adjust_prior`, and the defaults of
  `coloc.abf`.
+ Wallace C. Eliciting priors and relaxing the single causal variant assumption in colocalisation
  analyses. _PLoS Genetics_ 16(4):e1008720, 2020. Source of the $Q$ dependence of prior odds, the
  suggestion $p_(12) = 5 times 10^(-6)$, and the sensitivity reweighting.
+ `coloc_abf.py` on branch `large-scale-protein-coloc`.

*Cited but not read.*

#set enum(start: 5)
+ Wakefield J. Bayes factors for genome-wide association studies: comparison with P-values. _Genetic
  Epidemiology_ 33:79-86, 2009. The origin of the approximate Bayes factor, per [1]. It is closed access
  and was not obtained, so attributions to it rest on [1]. The derivation above does not rely on it.
+ Prentice RL, Pyke R. Logistic disease incidence models and case-control studies. _Biometrika_
  66:403-411, 1979. Cited for the classical result in @sec-logistic, stated from memory.

*What in this note is not in the sources.* The sources state or imply the results of @sec-abf and
@sec-posterior, though usually more briefly. @thm-abf's exactness for the Gaussian model and the
explanation of why LD is absent, @prop-implemented, @thm-reduction and its consequences in
@sec-structure, and @prop-masked were derived for this note and should be judged on their proofs.

*Numerical verification.* `verify_coloc_abf_derivation.py`, in the directory of this file, tests the
claims above against independent computation. @thm-abf agrees with multivariate normal densities under
LD and covariates to $2 times 10^(-13)$. @thm-posterior agrees with brute-force enumeration of all
$(Q+1)^2$ configurations using full joint densities, and that enumeration agrees with `ColocAbf` on
the branch to $10^(-14)$. @thm-reduction agrees with `combine_abf` on 2000 random inputs to
$10^(-14)$. @lem-prior agrees with enumeration of all $4^6$ category vectors. The script also
regenerates @tab-logistic, @tab-stylized, and the numerical values quoted in the text, apart from
simple arithmetic on them.
