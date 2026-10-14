# Why the SUSIE adjustment shifted on the numpy 1.26 -> 2.5 upgrade

Investigating the `test_susie_r_finemap_task` failure under the Python 3.13 upgrade,
where the null-model LD adjustment moved from `<= 0.01` to `0.0163` and broke the test's
sanity bound.

Reproduction: `experiments/claude/susie_numpy2_diag.py`, run in both envs. Logs in
`experiments/claude/logs/susie_numpy2_diag_{main_numpy1,branch_numpy2}.log`.

## Result: the R computation didn't change — its Python-computed inputs did

Fingerprints of each intermediate, same seed (`default_rng(40)`), causal_variants = []:

| quantity | numpy 1.26.4 | numpy 2.5.2 | |
| --- | --- | --- | --- |
| `raw_normals` (fresh PCG64 draw) | md5 `6b2b0689` | md5 `6b2b0689` | identical |
| `covar` | md5 `8dff0e80` | md5 `8dff0e80` | identical |
| `genotypes` (multivariate_normal) | md5 `6accf9aa` | md5 `ea3b79d1` | **differ** |
| `genotypes` sum-of-squares | `1.170593652119e+05` | `1.170593652119e+05` | **identical** |
| `corrcoef` / `beta` / `zscores` | ... | ... | differ (downstream) |
| adjustment (proxy) | 0.2717 | 0.2582 | differ |

The divergence originates precisely at `genotypes = generator.multivariate_normal(...)`:

- The raw bit-generator stream (`raw_normals`) is **identical** — numpy guarantees PCG64
  stream stability across versions, so it is NOT the RNG.
- The covariance matrix (`covar`) is identical.
- `genotypes` differ, **but their sum-of-squares is bit-identical.** That is the signature
  of a different-but-valid covariance factorization: `multivariate_normal` draws
  `genotypes = standard_normals @ A` where `A A^T = covar`. A different factor `A` gives a
  different sample, yet `sum((X@A)^2) = trace(X covar X^T)` is invariant to the choice of
  `A`. So the same random numbers were rotated by a different valid factor.

`multivariate_normal` factorizes `covar` with an **SVD** (the Generator default), computed
by LAPACK/OpenBLAS. The two numpy wheels bundle different OpenBLAS builds:

- numpy 1.26.4 wheel: **OpenBLAS 0.3.23**
- numpy 2.5.2 wheel: **OpenBLAS 0.3.34** (scipy-openblas)

Different OpenBLAS -> different (equally valid) SVD of `covar` -> a different genotype
realization -> different `corrcoef` -> different R `univariate_regression` betas/SEs ->
different z-scores -> different `estimate_s_rss` adjustment.

## Takeaways

1. This is **not a bug and not precision loss.** Both draws are correct samples from the
   same distribution; the pipeline still works (in the real test the causal-variant PIPs and
   credible sets are unaffected). Only the arbitrary `<= 0.01` sanity bound, calibrated to the
   old realization, is exceeded.
2. `Generator.multivariate_normal` is **not reproducible across numpy/BLAS versions** — only
   the underlying bit stream is guaranteed stable. The test implicitly relied on cross-version
   reproducibility that numpy never promised. A second BLAS-sensitive step,
   `make_psd_corr`'s `eigh`, would compound this if its `tol=1e-4` branch were near the edge
   (here it was not: min eigenvalue ~0.089, branch not taken in either env).
3. Options to make the test robust: (a) relax the bound (the substantive PIP/credible-set
   assertions still pass); or (b) make the synthetic data BLAS-independent — e.g. draw
   standard normals and multiply by a fixed Cholesky factor of `covar` computed once, rather
   than relying on `multivariate_normal`'s SVD factor. Cholesky is mathematically unique for a
   PD matrix, so it is far more stable across BLAS builds than SVD.

## Distribution of the adjustment (calibrating a realization-independent bound)

`experiments/claude/susie_adjustment_distribution.py` Monte-Carlos the adjustment over 300
independent draws of the null-model DGP, faithful to the real task (LD = full corrcoef via
M + M.T, then make_psd_corr). Seed 40 reproduces the real failing value exactly (0.016337).

    median = 0.000   mean = 0.027   p90 = 0.090   p95 = 0.117
    p99 = 0.181   p99.9 = 0.248   max = 0.265

The adjustment is NOT reliably small: the median is 0 but the right tail reaches ~0.27. The
old `<= 0.01` bound is violated by a large fraction of draws — it was never a valid
distribution-level bound, only an artifact of the one realization that happened to pass.

`estimate_s_rss` returns s in [0, 1), so a bound in [0.3, 0.4] still meaningfully guards
against a gross regression (s blowing up toward 1) while passing all legitimate draws:

- `<= 0.3` holds in >99.9% of the 300 draws (~1.13x the observed max) — recommended.
- `<= 0.2` holds in ~99% (p99 = 0.181) if a small documented failure rate is acceptable.

Cost/benefit: since the adjustment naturally ranges 0–0.27, this assertion is a weak guard;
the substantive assertions are PIP >= 0.95 and the credible-set count. The adjustment bound
is worth keeping only as a loose blow-up guard (<= 0.3), or dropping in favor of the
signal-based assertions.

## Why the adjustment is large even though z and R come from the same data

`experiments/claude/susie_adjustment_mechanism.py`. estimate_s_rss (null-mle) fits
z ~ N(0, (1-s)R + sI) and returns s in [0,1). In R's eigenbasis (R = V diag(d) V^T,
u = V^T z), correct specification gives u_k ~ N(0, d_k), so the per-direction "leverage"
u_k^2 / d_k is ~chi-square(1): mean 1 but variance ~ ... dominated by the small-d_k
directions. The MLE inflates s whenever z carries above-expected leverage in R's
small-eigenvalue directions, because the likelihood penalises under-modelled variance there
most heavily.

This LD matrix is ill-conditioned (the SNPs share a strong common factor q, cond(R) ~ 260),
so the small-eigenvalue directions have high-variance leverage. Evidence over 200 draws:

    frac(s == 0) varies by seed range; mean_s ~ 0.03, max_s ~ 0.21
    Spearman(s, summed leverage in 10 smallest-eigenvalue dirs) = +0.57   <-- strong
    Spearman(s, max_k u_k^2/d_k)                                = +0.31
    Spearman(s, min eigenvalue of R)                            = -0.04   <-- ~none

Interpretation: s is NOT detecting a systematic z/R mismatch — there is none (median s ~ 0,
no bias direction; z ~ N(0,R) holds in distribution because both come from the same data).
s is a one-sided, noise-driven diagnostic of how atypical THIS realization of z is in R's
smallest-eigenvalue directions. Because R is ill-conditioned, those directions fluctuate a
lot, so s occasionally lands well above 0. The near-zero correlation with min(eig) (the
spectrum is ~fixed across draws) vs the strong correlation with the random per-draw leverage
confirms it is the projection of z, not the conditioning per se, that drives each draw.
