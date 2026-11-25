#set document(title: "How _fwl_target_coef applies the Frisch–Waugh–Lovell theorem")
#set page(paper: "us-letter", margin: 1in, numbering: "1")
#set text(size: 11pt)
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")
#show raw: set text(size: 9.5pt)

#align(center)[
  #text(size: 16pt, weight: "bold")[How `_fwl_target_coef` applies the Frisch–Waugh–Lovell theorem]
  #v(0.3em)
  Companion note to the vectorized S-LDSC plan, Task 5 (shared_controls mode)
]

= What the function computes

For one control group, the shared_controls mode fits $k$ separate weighted least squares (WLS) regressions, one per target annotation $j = 1, dots, k$:

$ y = x_j beta_j + Z delta_j + e_j, quad "weighted by" W = "diag"(w). $ <model>

The pieces are as follows.

- $y in RR^S$ is the vector of $chi^2$ statistics over the $S$ regression SNPs (`inputs.chi2`).
- $Z in RR^(S times q)$ is the shared design (`shared.design`, called `a` / $A$ in the code). Its $q$ columns are the intercept, the baseline LD-score columns and the control-group LD-score column, each scaled by $N slash overline(N)$.
- $x_j in RR^S$ is the scaled LD score of target $j$ (column $j$ of `t` in the code).
- $W = "diag"(w)$ is the $S times S$ diagonal weight matrix. In shared_controls mode, $w$ is computed once per control group from the model _without_ a target column, so the same $W$ is used for all $k$ targets. This is the only approximation in the mode. Given $W$, everything below is exact algebra.
- $beta_j$ is the scalar target coefficient. It is the only quantity we want; $delta_j in RR^q$ is a nuisance.

`_fwl_target_coef` returns $hat(beta)_1, dots, hat(beta)_k$. It does so without ever forming the $(q+1) times (q+1)$ normal equations of each regression. It only uses inner products that are shared across the $k$ targets or are cheap per target.

= The theorem as stated on Wikipedia

Wikipedia states the theorem for the ordinary least squares (OLS) regression

$ y = X hat(beta) + Z hat(delta) + hat(e), $

where $X$ is the $n times k$ matrix of regressors of interest and $Z$ is the $n times p$ matrix of other regressors. (Wikipedia writes transposes as $X'$; this note writes $X^top$.) Define the projection and annihilator matrices

$ P_Z = Z (Z^top Z)^(-1) Z^top, quad M_Z = I - P_Z, $

and the residualized quantities

$ tilde(X) = M_Z X = X - Z (Z^top Z)^(-1) Z^top X, quad tilde(y) = M_Z y = y - Z (Z^top Z)^(-1) Z^top y. $

The theorem says that the coefficient on $X$ in the full regression equals the coefficient from the "double residual" regression of $tilde(y)$ on $tilde(X)$:

$ hat(beta) = (tilde(X)^top tilde(X))^(-1) tilde(X)^top tilde(y). $ <fwl>

In words, partial $Z$ out of both $y$ and $X$, then regress residual on residual.

= From OLS to WLS

WLS with weights $W$ is OLS on rescaled data. Let $W^(1 slash 2) = "diag"(sqrt(w))$ and define

$ y^* = W^(1 slash 2) y, quad x_j^* = W^(1 slash 2) x_j, quad Z^* = W^(1 slash 2) Z. $

Minimizing $(y - x_j beta - Z delta)^top W (y - x_j beta - Z delta)$ is the same as minimizing $norm(y^* - x_j^* beta - Z^* delta)^2$. So we can apply @fwl to the starred data, with $X = x_j^*$ (a single column, so Wikipedia's $k = 1$) and with $Z^*$ in the role of $Z$. Every inner product of starred quantities is a $W$-weighted inner product of the originals, for example

$ Z^*^top Z^* = Z^top W Z, quad Z^*^top y^* = Z^top W y, quad x_j^*^top x_j^* = x_j^top W x_j. $

= Expanding the residual inner products

Since $x_j^*$ is one column, $tilde(X)^top tilde(X)$ is a scalar and @fwl becomes a ratio. Write $M^* = I - Z^* (Z^*^top Z^*)^(-1) Z^*^top$. $M^*$ is symmetric and idempotent ($M^*^top M^* = M^*$), so

$ hat(beta)_j = (tilde(x)_j^top tilde(y)) / (tilde(x)_j^top tilde(x)_j) = (x_j^*^top M^* y^*) / (x_j^*^top M^* x_j^*). $

Expanding $M^*$ and converting the starred inner products back to weighted ones gives

$ hat(beta)_j = (x_j^top W y - x_j^top W Z (Z^top W Z)^(-1) Z^top W y) / (x_j^top W x_j - x_j^top W Z (Z^top W Z)^(-1) Z^top W x_j). $ <ratio>

This is the formula in the docstring of `_fwl_target_coef`, written there with $A$ for $Z$ and $t_c$ for $x_j$. Every term in @ratio is built from five inner products, which are exactly the function's five arguments:

#table(
  columns: 4,
  align: (left, left, left, left),
  stroke: 0.5pt,
  [*argument*], [*quantity*], [*shape (full data)*], [*shared across targets?*],
  [`a_w_a`], [$Z^top W Z$], [$(q, q)$], [yes],
  [`a_w_y`], [$Z^top W y$], [$(q)$], [yes],
  [`a_w_t`], [$Z^top W x_j$, stacked over $j$], [$(q, k)$], [no],
  [`t_w_t`], [$x_j^top W x_j$, stacked over $j$], [$(k)$], [no],
  [`t_w_y`], [$x_j^top W y$, stacked over $j$], [$(k)$], [no],
)

Collect the targets into $X = [x_1, dots, x_k] in RR^(S times k)$ (the chunk `t`), and write $G = Z^top W Z$, $g = Z^top W y$ and $C = Z^top W X$ (so column $j$ of $C$ is $c_j = Z^top W x_j$).

= Line-by-line correspondence with the code

```python
solution = np.linalg.solve(a_w_a, np.concatenate([a_w_y[..., None], a_w_t], axis=-1))
```

This solves $G thin [s_y, s_1, dots, s_k] = [g, c_1, dots, c_k]$ in a single call. One factorization of $G$ serves all $k + 1$ right-hand sides:

$ s_y = G^(-1) g = (Z^top W Z)^(-1) Z^top W y, quad s_j = G^(-1) c_j = (Z^top W Z)^(-1) Z^top W x_j. $

These are the coefficients of the auxiliary regressions of $y$ on $Z$ and of $x_j$ on $Z$. They are the "residualize on $Z$" step of the theorem. `solution[..., :1]` is $s_y$ and `solution[..., 1:]` is $S_X = [s_1, dots, s_k]$, of shape $(q, k)$.

```python
numerator = t_w_y - (a_w_t * solution[..., :1]).sum(axis=-2)
```

`a_w_t * solution[..., :1]` multiplies each column $c_j$ elementwise by $s_y$ (broadcast over $k$). Summing over the $q$ axis (`axis=-2`) gives $c_j^top s_y$. Because $G$ is symmetric,

$ c_j^top s_y = x_j^top W Z (Z^top W Z)^(-1) Z^top W y, $

so `numerator[j]` is the numerator of @ratio:

$ "numerator"_j = x_j^top W y - x_j^top W Z (Z^top W Z)^(-1) Z^top W y = tilde(x)_j^top tilde(y). $

```python
denominator = t_w_t - (a_w_t * solution[..., 1:]).sum(axis=-2)
```

`a_w_t * solution[..., 1:]` is the elementwise product of $C$ and $S_X$, both $(q, k)$. Summing over $q$ gives $c_j^top s_j$ for each $j$, which is

$ c_j^top s_j = x_j^top W Z (Z^top W Z)^(-1) Z^top W x_j. $

This is the diagonal of $C^top G^(-1) C$, computed without forming the off-diagonal entries. Those entries would couple different targets, but each target sits in its own regression, so they are never needed. Therefore

$ "denominator"_j = x_j^top W x_j - x_j^top W Z (Z^top W Z)^(-1) Z^top W x_j = tilde(x)_j^top tilde(x)_j. $

```python
return numerator / denominator
```

This returns $hat(beta)_j = tilde(x)_j^top tilde(y) slash tilde(x)_j^top tilde(x)_j$ for every $j$, which is @fwl for a single regressor of interest.

= Why this is cheaper than solving each regression

Solving @model directly needs the $(q + 1) times (q + 1)$ Gram matrix $[x_j, Z]^top W [x_j, Z]$ and one solve per target. In exact mode, $W$ also differs per target, which forces a fresh $Z^top W_j Z$ for each one. With a shared $W$:

- $Z^top W Z$ and $Z^top W y$ are formed once per control group (and per jackknife block). They do not depend on $j$.
- Each target only adds $c_j = Z^top W x_j$ ($q$ numbers) and the two scalars $x_j^top W x_j$ and $x_j^top W y$.
- One factorization of the $q times q$ matrix $G$ handles all $k$ targets.

The theorem is what lets us drop the $(q+1)$-sized system per target: the target coefficient depends on $Z$ only through the shared quantities $G$ and $g$ and the per-target vector $c_j$.

= Jackknife delete values: broadcasting over blocks

`_shared_controls_group` computes every inner product block by block. SNPs are split into contiguous jackknife blocks $b = 1, dots, n_"blocks"$ by the separators. Because $W$ is diagonal, each inner product is a sum over SNPs, so it is a sum over blocks:

$ Z^top W Z = sum_(b) Z_b^top W_b Z_b, quad "and likewise for" Z^top W y, thin Z^top W x_j, thin x_j^top W x_j, thin x_j^top W y. $

The per-block terms are stacked along a leading axis of length $n_"blocks"$: `a_w_a` has shape $(n_"blocks", q, q)$, and so on. Then:

- *Full-data estimate.* Passing the sums over blocks (`.sum(axis=0)`) gives arguments with no leading axis, and @ratio yields $hat(beta)_j$ with shape $(k)$.
- *Delete-block-$b$ estimate.* Leaving block $b$ out of the regression removes exactly its term from every sum, so the leave-one-out inner products are "total minus block $b$". This is `x.sum(axis=0)[None] - x`, with shape $(n_"blocks", ...)$. `np.linalg.solve` and the elementwise products broadcast over that leading axis. One call therefore returns all delete values $hat(beta)_j^((-b))$ with shape $(n_"blocks", k)$, each the FWL ratio @ratio evaluated on the data with block $b$ removed.

The weights $w$ are not recomputed when a block is deleted. This matches ldsc's block jackknife, which deletes rows of the already-weighted design.

= Scope of the approximation

The FWL step introduces no error: for a given $W$, @ratio is algebraically identical to the target coefficient of the full $(q + 1)$-column WLS fit. The plan's dense-WLS oracle test checks this at rtol $10^(-8)$. The difference between shared_controls and exact mode comes entirely from the choice of $W$:

- *exact* recomputes the ldsc heteroskedasticity weights for each target model. Those weights depend on $h^2_"agg"$ and on the total LD score, both of which include the target column.
- *shared_controls* reuses the weights of the target-free model.

That is the accuracy cost that Task 7's real-data comparison measures.
