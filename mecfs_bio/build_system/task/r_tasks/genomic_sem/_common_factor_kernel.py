"""
Vectorized common-factor multivariate GWAS.

This module fits, jointly across many SNPs, a structural equation model in
which a single latent factor F1 explains genetic covariance among k = 2
phenotypes (traits) and is itself regressed on a single SNP. The quantity of
scientific interest is the SNP -> F1 regression coefficient and its standard
error, computed once per SNP across the genome.

It is a numpy/scipy reimplementation of GenomicSEM's `commonfactorGWAS`
routine; we replace its per-SNP `lavaan::sem` call with a batched
Gauss-Newton DWLS optimization, eliminating ~20 ms of lavaan dispatch
overhead per SNP. v1 is restricted to k = 2 traits and `GC = "standard"`;
the math generalizes cleanly.

References:
- Grotzinger et al. 2019, Nat Hum Behav 3:513-525  (GenomicSEM paper)
- GenomicSEM wiki, page 4 ("Common Factor GWAS")
- GenomicSEM::commonfactorGWAS R source, inspected via deparse()


# Background and inputs

For each of the k traits, a single-trait GWAS has previously been run,
producing per-SNP marginal regression coefficients `beta_SNP[i, j]` and
standard errors `SE_SNP[i, j]` (i indexes SNPs, j indexes traits).

LD score regression (LDSC, Bulik-Sullivan et al. 2015) summarizes those
sumstats into three trait-level matrices:

- `S_LD`  (k, k):  the *genetic* covariance matrix among traits.
                   `S_LD[j, j]` is heritability of trait j; off-diagonals are
                   genetic covariances.
- `V_LD`  (m, m), m = k(k+1)/2:  the sampling covariance of `vech(S_LD)` -
                   i.e. the second-moment uncertainty in S_LD's entries
                   given the finite GWAS sample sizes.
- `I_LD`  (k, k):  LDSC intercepts. The diagonal absorbs uncontrolled
                   confounding; off-diagonals capture cross-trait sample
                   overlap (between-trait sampling correlation in the GWAS).

For each SNP we additionally need

- `varSNP[i] = 2 * MAF[i] * (1 - MAF[i])`  -- the genotype variance under
  Hardy-Weinberg equilibrium, where MAF is the minor allele frequency.
  Throughout, "SNP-trait covariance" means `varSNP * beta_SNP` (covariance
  on the genotype-dosage scale, given that GWAS betas are slopes of the
  trait on the standardized genotype).


# The structural equation model

We fit, separately for each SNP, the model

    F1 =~ trait_1 + lambda * trait_2             (factor loadings)
    F1  ~ beta * SNP                             (latent regression)
    trait_j ~ 0 * SNP                            (no direct SNP effect)

with free parameters `theta = (lambda, beta, psi_1, psi_2, sigma_sq)`:

- lambda   :  factor loading on trait 2 (loading on trait 1 fixed to 1 for
              identification; this is lavaan's default scaling).
- beta     :  the SNP -> F1 regression coefficient.  *This is the quantity
              we report.*
- psi_1,   :  residual variances of traits 1 and 2 (the part of each trait
  psi_2       not mediated by F1).
- sigma_sq :  residual variance of F1 after partialing out the SNP.

The model implies an augmented 3x3 covariance matrix over the variables
(SNP, trait_1, trait_2). Define `varF = beta^2 * varSNP + sigma_sq`, the
unconditional variance of F1. Stacking the lower triangle column-major
(lavaan's `vech` convention), the model-implied vector
`sigma(theta) in R^6` has entries

    sigma_1 = varSNP                                (Var(SNP),   constant)
    sigma_2 = beta * varSNP                         (Cov(SNP, t1))
    sigma_3 = lambda * beta * varSNP                (Cov(SNP, t2))
    sigma_4 = beta^2 * varSNP + sigma_sq + psi_1    (Var(t1))
    sigma_5 = lambda * varF                         (Cov(t1, t2))
    sigma_6 = lambda^2 * varF + psi_2               (Var(t2))

The corresponding *observed* vector `s in R^6` is `vech(S_Fullrun)`, where
`S_Fullrun` augments `S_LD` with a SNP row/column built from `varSNP` and
`beta_SNP`:

    s_1 = varSNP                                    (constant per-SNP)
    s_2 = varSNP * beta_SNP[i, 1]
    s_3 = varSNP * beta_SNP[i, 2]
    s_4 = S_LD[1, 1]                                (constant across SNPs)
    s_5 = S_LD[2, 1]                                (constant)
    s_6 = S_LD[2, 2]                                (constant)

Note that only entries 2 and 3 of `s` carry per-SNP information; the rest
are constants of the trait pair. The point of the per-SNP fit is to find
the `theta` that best reconciles those two SNP-trait covariances with the
trait-trait genetic covariance structure.


# DWLS objective and per-SNP weights

Both `s` and `sigma` are subject to sampling noise inherited from the
underlying GWAS sumstats. GenomicSEM's diagonally-weighted least squares
(DWLS) loss weights the squared residual on each vech entry by the inverse
of that entry's *marginal* sampling variance:

    L(theta) = sum_p (1 / V_Full[p, p]) * (s_p - sigma_p(theta))^2

where `V_Full` is the (k+1)(k+2)/2 = 6-dimensional sampling covariance of
`vech(S_Fullrun)`. It has block-diagonal structure (off-diagonal blocks
zero):

    V_Full = blockdiag(varSNPSE2,  V_SNP(SE_SNP, I_LD, varSNP),  V_LD)

- `varSNPSE2`  -- the (assumed-known) variance of the SNP variance estimate.
                  GenomicSEM uses `(5e-4)^2` by default.
- `V_SNP`  (k, k)  -- per-SNP sampling covariance of the SNP-trait
                       covariance row, derived from `SE_SNP` and `I_LD`.
                       Diagonal: `(SE_SNP[j] * sqrt(I_LD[j, j]) * varSNP)^2`.
                       Off-diagonal includes the cross-trait LDSC intercept
                       to account for sample overlap.
- `V_LD`  (m, m)   -- the LDSC-supplied sampling covariance of the
                       trait-trait block (does not vary per SNP).

`V_Full` *enters the optimization only via its diagonal*, since the loss is
DWLS. Its full matrix structure is needed only for the sandwich variance
below.


# Algorithm

The Gauss-Newton step for L(theta) is, with `J = d sigma / d theta` (a 6x5
Jacobian computed in closed form in `model_jacobian`) and weight vector
`w = diag(V_Full)^-1`:

    delta = (J' diag(w) J)^-1  J' diag(w) (s - sigma(theta))
    theta_new = theta + delta

Because the model structure is identical across SNPs and only the data
arrays (`varSNP`, `beta_SNP`, `SE_SNP`) vary per SNP, the entire iteration
broadcasts cleanly across a leading "SNP" axis. We start from
moment-of-method values (lambda_0 = S_LD[2,1] / S_LD[1,1], variances split
as in `_starting_theta`) and stop when `||delta||_inf < tol`. In practice
~5 iterations suffice.


# Sandwich variance

DWLS with weights `diag(V)^-1` is misspecified relative to the optimal
full-V GLS, so a sandwich correction is needed for valid inference. With
`A = J' diag(w) J` (the bread inverse) and `B = J' diag(w) V_Full diag(w) J`
(the meat),

    Var(theta_hat) ~ A^-1 B A^-1

and we report `se_c = sqrt(Var(theta_hat)[beta, beta])` along with
`Z = beta_hat / se_c` and a two-sided normal p-value.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from scipy.stats import norm

# Number of traits this v1 supports. Asserted on entry; generalising is
# straightforward but kept out of v1 to keep the math compact.
_K = 2

# Indices into theta = (lambda, beta, psi_1, psi_2, sigma_sq).
_LAMBDA_IDX = 0
_BETA_IDX = 1
_PSI_1_IDX = 2
_PSI_2_IDX = 3
_SIGMA_SQ_IDX = 4
_N_PARAMS = 5

# Number of vech entries for k+1 = 3 augmented variables: 3*4/2 = 6.
_VECH_LEN = 6

# Floor for variance parameters during Newton iterations; prevents zero/negative
# values from breaking the model-implied covariance computation.
_VAR_FLOOR = 1e-12


class CommonFactorGWASResult(NamedTuple):
    """Per-SNP results from `fit_common_factor_gwas`."""

    est: np.ndarray  # (N,) point estimate of beta
    se_c: np.ndarray  # (N,) sandwich SE for beta
    z: np.ndarray  # (N,) Z-statistic
    p: np.ndarray  # (N,) two-sided p-value
    fail: np.ndarray  # (N,) 1 if optimization or SE computation failed
    warning: np.ndarray  # (N,) 1 if non-PD V_Full / S_Fullrun (still attempted)
    n_iter: np.ndarray  # (N,) Newton iterations used


def compute_v_snp_batch(
    SE_SNP: np.ndarray, I_LD: np.ndarray, varSNP: np.ndarray
) -> np.ndarray:
    """
    Vectorized GenomicSEM::.get_V_SNP for GC="standard".

    SE_SNP   : (N, k)
    I_LD     : (k, k)
    varSNP   : (N,)

    Returns (N, k, k).
    """
    N, k = SE_SNP.shape
    assert k == _K, f"_common_factor_kernel only supports k=2 (got k={k})"
    assert I_LD.shape == (k, k)
    assert varSNP.shape == (N,)

    sqrt_diag_I = np.sqrt(np.diag(I_LD))  # (k,)
    # Scale per-trait SE by sqrt(diagonal LDSC intercept)
    se_scaled = SE_SNP * sqrt_diag_I  # (N, k)
    # Outer product per SNP: (N, k, k)
    se_outer = se_scaled[:, :, None] * se_scaled[:, None, :]
    # Element-wise factor: I_LD[x,y] for off-diagonal, 1 for diagonal so the
    # diagonal works out to (SE * sqrt(I)) * (SE * sqrt(I)) * varSNP^2.
    intercept_factor = I_LD.copy()
    np.fill_diagonal(intercept_factor, 1.0)
    var_snp_sq = varSNP**2  # (N,)
    return se_outer * intercept_factor[None, :, :] * var_snp_sq[:, None, None]


def compute_v_full_batch(
    V_LD: np.ndarray, varSNPSE2: float, V_SNP_batch: np.ndarray
) -> np.ndarray:
    """
    Vectorized GenomicSEM::.get_V_full.

    V_LD         : (3, 3) for k=2 (sampling covariance of vech(S_LD))
    varSNPSE2    : scalar (default (5e-4)^2)
    V_SNP_batch  : (N, k, k) from compute_v_snp_batch

    Returns (N, M, M) where M = (k+1)(k+2)/2 = 6.
    """
    N = V_SNP_batch.shape[0]
    k = V_SNP_batch.shape[1]
    assert k == _K
    M = (k + 1) * (k + 2) // 2  # 6
    assert V_LD.shape == (k * (k + 1) // 2, k * (k + 1) // 2)

    V_full = np.zeros((N, M, M))
    V_full[:, 0, 0] = varSNPSE2
    V_full[:, 1 : 1 + k, 1 : 1 + k] = V_SNP_batch
    V_full[:, 1 + k : M, 1 + k : M] = V_LD[None, :, :]
    return V_full


def compute_s_fullrun_batch(
    S_LD: np.ndarray, varSNP: np.ndarray, beta_SNP: np.ndarray
) -> np.ndarray:
    """
    Vectorized augmented sample covariance over (SNP, trait_1, ..., trait_k).

    S_LD     : (k, k)
    varSNP   : (N,)
    beta_SNP : (N, k)

    Returns (N, k+1, k+1).
    """
    N, k = beta_SNP.shape
    assert k == _K
    assert S_LD.shape == (k, k)
    assert varSNP.shape == (N,)

    S = np.zeros((N, k + 1, k + 1))
    S[:, 0, 0] = varSNP
    snp_trait_cov = varSNP[:, None] * beta_SNP  # (N, k)
    S[:, 1:, 0] = snp_trait_cov
    S[:, 0, 1:] = snp_trait_cov
    S[:, 1:, 1:] = S_LD[None, :, :]
    return S


def vech_lower_col_major(M: np.ndarray) -> np.ndarray:
    """
    vech of lower triangle in column-major order (lavaan convention).

    M : (..., n, n)
    Returns (..., n*(n+1)/2)
    """
    n = M.shape[-1]
    rows: list[int] = []
    cols: list[int] = []
    for j in range(n):
        for i in range(j, n):
            rows.append(i)
            cols.append(j)
    return M[..., rows, cols]


def model_implied_vech(theta: np.ndarray, varSNP: np.ndarray) -> np.ndarray:
    """
    Compute vech(Sigma(theta)) per SNP for the common-factor model.

    theta   : (N, 5) columns (lambda, beta, psi_1, psi_2, sigma_sq)
    varSNP  : (N,)

    Returns (N, 6).
    """
    lam = theta[:, _LAMBDA_IDX]
    beta = theta[:, _BETA_IDX]
    psi_1 = theta[:, _PSI_1_IDX]
    psi_2 = theta[:, _PSI_2_IDX]
    sigma_sq = theta[:, _SIGMA_SQ_IDX]
    var_F = beta**2 * varSNP + sigma_sq

    sigma = np.empty((theta.shape[0], _VECH_LEN))
    sigma[:, 0] = varSNP
    sigma[:, 1] = beta * varSNP
    sigma[:, 2] = lam * beta * varSNP
    sigma[:, 3] = beta**2 * varSNP + sigma_sq + psi_1
    sigma[:, 4] = lam * var_F
    sigma[:, 5] = lam**2 * var_F + psi_2
    return sigma


def model_jacobian(theta: np.ndarray, varSNP: np.ndarray) -> np.ndarray:
    """
    Compute J = d vech(Sigma) / d theta per SNP. Shape (N, 6, 5).

    Rows correspond to vech entries (SNP-SNP, SNP-t1, SNP-t2, t1-t1, t1-t2, t2-t2).
    Columns correspond to (lambda, beta, psi_1, psi_2, sigma_sq).
    """
    N = theta.shape[0]
    lam = theta[:, _LAMBDA_IDX]
    beta = theta[:, _BETA_IDX]
    sigma_sq = theta[:, _SIGMA_SQ_IDX]
    v = varSNP
    var_F = beta**2 * v + sigma_sq

    J = np.zeros((N, _VECH_LEN, _N_PARAMS))

    # Row 0 (sigma_1 = varSNP) -> all zero, no parameters affect it.

    # Row 1 (sigma_2 = beta * varSNP)
    J[:, 1, _BETA_IDX] = v

    # Row 2 (sigma_3 = lambda * beta * varSNP)
    J[:, 2, _LAMBDA_IDX] = beta * v
    J[:, 2, _BETA_IDX] = lam * v

    # Row 3 (sigma_4 = beta^2 * varSNP + sigma_sq + psi_1)
    J[:, 3, _BETA_IDX] = 2.0 * beta * v
    J[:, 3, _PSI_1_IDX] = 1.0
    J[:, 3, _SIGMA_SQ_IDX] = 1.0

    # Row 4 (sigma_5 = lambda * var_F)
    J[:, 4, _LAMBDA_IDX] = var_F
    J[:, 4, _BETA_IDX] = 2.0 * lam * beta * v
    J[:, 4, _SIGMA_SQ_IDX] = lam

    # Row 5 (sigma_6 = lambda^2 * var_F + psi_2)
    J[:, 5, _LAMBDA_IDX] = 2.0 * lam * var_F
    J[:, 5, _BETA_IDX] = 2.0 * lam**2 * beta * v
    J[:, 5, _PSI_2_IDX] = 1.0
    J[:, 5, _SIGMA_SQ_IDX] = lam**2

    return J


def _starting_theta(
    S_LD: np.ndarray, varSNP: np.ndarray, beta_SNP: np.ndarray
) -> np.ndarray:
    """
    Per-SNP starting values for Newton iteration.

    Returns (N, 5).
    """
    N = varSNP.shape[0]
    s11 = S_LD[0, 0]
    s21 = S_LD[1, 0]
    s22 = S_LD[1, 1]

    lam0 = s21 / s11 if s11 > 0 else 1.0
    sigma_sq_0 = max(_VAR_FLOOR, s11 * 0.5)
    psi_1_0 = max(_VAR_FLOOR, s11 * 0.5)
    psi_2_0 = max(_VAR_FLOOR, s22 - lam0 * lam0 * sigma_sq_0)

    theta = np.empty((N, _N_PARAMS))
    theta[:, _LAMBDA_IDX] = lam0
    theta[:, _BETA_IDX] = beta_SNP[:, 0]  # OLS estimate of trait 1 regression
    theta[:, _PSI_1_IDX] = psi_1_0
    theta[:, _PSI_2_IDX] = psi_2_0
    theta[:, _SIGMA_SQ_IDX] = sigma_sq_0
    return theta


def _project_to_feasible(theta: np.ndarray) -> None:
    """In-place project to feasible region (variances >= floor)."""
    np.maximum(theta[:, _PSI_1_IDX], _VAR_FLOOR, out=theta[:, _PSI_1_IDX])
    np.maximum(theta[:, _PSI_2_IDX], _VAR_FLOOR, out=theta[:, _PSI_2_IDX])
    np.maximum(theta[:, _SIGMA_SQ_IDX], _VAR_FLOOR, out=theta[:, _SIGMA_SQ_IDX])


def _gauss_newton_fit(
    S_LD: np.ndarray,
    V_full: np.ndarray,
    varSNP: np.ndarray,
    beta_SNP: np.ndarray,
    *,
    max_iter: int,
    tol: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Vectorized Gauss-Newton DWLS optimization for the common-factor model.

    Returns (theta, n_iter, converged) where:
    - theta     : (N, 5) final parameter estimates
    - n_iter    : (N,) Newton iterations used
    - converged : (N,) bool flag
    """
    N = varSNP.shape[0]

    # Augmented sample vech(S_Fullrun)
    S_full = compute_s_fullrun_batch(S_LD, varSNP, beta_SNP)
    s_obs = vech_lower_col_major(S_full)  # (N, 6)

    # Diagonal weights: 1 / diag(V_full) per SNP. Note: V_full[..., 0, 0]
    # corresponds to sigma_1 which has zero Jacobian, so its weight value
    # has no influence on the fit -- but it matters for the sandwich SE.
    W_diag = 1.0 / np.diagonal(V_full, axis1=-2, axis2=-1)  # (N, 6)

    theta = _starting_theta(S_LD, varSNP, beta_SNP)
    _project_to_feasible(theta)

    n_iter = np.zeros(N, dtype=np.int32)
    converged = np.zeros(N, dtype=bool)
    active = np.ones(N, dtype=bool)

    for it in range(max_iter):
        if not active.any():
            break
        idx = np.where(active)[0]
        sub_theta = theta[idx]
        sub_var = varSNP[idx]
        sub_obs = s_obs[idx]
        sub_W = W_diag[idx]

        sigma = model_implied_vech(sub_theta, sub_var)
        residual = sub_obs - sigma  # (n_active, 6)
        J = model_jacobian(sub_theta, sub_var)  # (n_active, 6, 5)

        # Gauss-Newton normal equations: (J' diag(W) J) delta = J' diag(W) r
        WJ = sub_W[:, :, None] * J  # (n_active, 6, 5)
        H = np.einsum("npj,npk->njk", J, WJ)  # (n_active, 5, 5)
        rhs = np.einsum("npj,np->nj", WJ, residual)  # (n_active, 5)

        # Solve per-SNP. np.linalg.solve raises if any matrix is singular --
        # catch and mark per-SNP via a manual solve loop with a try/except is
        # too slow. Instead, use np.linalg.solve with regularisation: try, and
        # if it fails, fall back to least-squares per SNP.
        try:
            delta = np.linalg.solve(H, rhs)  # (n_active, 5)
        except np.linalg.LinAlgError:
            # Fall back to per-SNP solve so we can isolate failures
            delta = np.zeros_like(rhs)
            for j in range(idx.size):
                try:
                    delta[j] = np.linalg.solve(H[j], rhs[j])
                except np.linalg.LinAlgError:
                    delta[j] = np.nan

        new_theta = sub_theta + delta
        # Project to feasible region (clip variances)
        _project_to_feasible(new_theta)

        # Detect divergence
        bad = ~np.all(np.isfinite(new_theta), axis=1)
        new_theta[bad] = sub_theta[bad]  # keep last value to bail out

        # Convergence test on the (clipped) update
        achieved_delta = new_theta - sub_theta
        max_step = np.max(np.abs(achieved_delta), axis=1)
        sub_converged = (max_step < tol) & ~bad

        theta[idx] = new_theta
        n_iter[idx] = it + 1
        converged[idx[sub_converged]] = True
        # Mark diverged SNPs as no longer active and not converged
        new_active = active.copy()
        new_active[idx[sub_converged | bad]] = False
        active = new_active

    return theta, n_iter, converged


def _sandwich_se(
    theta: np.ndarray, varSNP: np.ndarray, V_full: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sandwich SE for theta given fitted estimates.

    Returns (se_per_param, fail_mask):
    - se_per_param: (N, 5) standard errors for each parameter
    - fail_mask: (N,) bool indicating which SNPs had numerical failure
    """
    N = theta.shape[0]
    J = model_jacobian(theta, varSNP)  # (N, 6, 5)
    W_diag = 1.0 / np.diagonal(V_full, axis1=-2, axis2=-1)  # (N, 6)

    WJ = W_diag[:, :, None] * J
    A = np.einsum("npj,npk->njk", J, WJ)  # bread^-1 = J' W J  (N, 5, 5)

    # Meat = J' W V_full W J
    WVW = W_diag[:, :, None] * V_full * W_diag[:, None, :]  # (N, 6, 6)
    meat = np.einsum("npj,npq,nqk->njk", J, WVW, J)  # (N, 5, 5)

    fail = np.zeros(N, dtype=bool)
    var_theta = np.full((N, _N_PARAMS, _N_PARAMS), np.nan)
    for n in range(N):
        try:
            bread = np.linalg.inv(A[n])
        except np.linalg.LinAlgError:
            fail[n] = True
            continue
        var_theta[n] = bread @ meat[n] @ bread

    diag_var = np.diagonal(var_theta, axis1=-2, axis2=-1).copy()  # (N, 5)
    bad = ~np.all(np.isfinite(diag_var), axis=1) | np.any(diag_var < 0, axis=1)
    fail = fail | bad

    se = np.sqrt(np.where(diag_var < 0, np.nan, diag_var))
    return se, fail


def fit_common_factor_gwas(
    *,
    S_LD: np.ndarray,
    V_LD: np.ndarray,
    I_LD: np.ndarray,
    beta_SNP: np.ndarray,
    SE_SNP: np.ndarray,
    varSNP: np.ndarray,
    varSNPSE2: float = (5e-4) ** 2,
    max_iter: int = 25,
    tol: float = 1e-10,
) -> CommonFactorGWASResult:
    """
    Per-SNP common-factor GWAS via vectorized DWLS Gauss-Newton + sandwich SE.

    Inputs (k = number of traits, currently must be 2; N = number of SNPs):
        S_LD       : (k, k)   trait-trait genetic covariance from LDSC
        V_LD       : (m, m) where m = k(k+1)/2 = 3 for k=2.  Sampling
                            covariance of vech(S_LD) from LDSC.
        I_LD       : (k, k)   LDSC intercepts (cross-trait sampling correlation)
        beta_SNP   : (N, k)   per-SNP standardised effect sizes per trait
        SE_SNP     : (N, k)   per-SNP standard errors per trait
        varSNP     : (N,)     2 * MAF * (1 - MAF) per SNP
        varSNPSE2  : scalar   variance of SNP-SNP entry sampling error
                              (default (5e-4)^2 matching GenomicSEM's SNPSE)
        max_iter, tol         Gauss-Newton convergence parameters

    Returns CommonFactorGWASResult with per-SNP arrays of length N.
    """
    N, k = beta_SNP.shape
    assert k == _K, f"v1 only supports k=2 traits; got k={k}"
    assert SE_SNP.shape == (N, k)
    assert varSNP.shape == (N,)
    assert S_LD.shape == (k, k)
    assert I_LD.shape == (k, k)
    m = k * (k + 1) // 2
    assert V_LD.shape == (m, m), f"V_LD must be ({m},{m}); got {V_LD.shape}"

    # Build V_SNP_batch and V_full per SNP (matches GenomicSEM internals).
    V_SNP_batch = compute_v_snp_batch(SE_SNP, I_LD, varSNP)
    V_full = compute_v_full_batch(V_LD, varSNPSE2, V_SNP_batch)

    # PD check for V_full (purely diagnostic; Gauss-Newton uses only the
    # diagonal of V_full for optimization). Mark non-PD SNPs with warning=1.
    eig_min_V = np.linalg.eigvalsh(V_full).min(axis=-1)
    warn_v = eig_min_V <= 0

    # Optimization
    theta, n_iter, converged = _gauss_newton_fit(
        S_LD, V_full, varSNP, beta_SNP, max_iter=max_iter, tol=tol
    )

    # Sandwich SE
    se_all, se_fail = _sandwich_se(theta, varSNP, V_full)

    fail = (~converged) | se_fail | (~np.all(np.isfinite(theta), axis=1))
    warning = warn_v.astype(np.int8)

    est = theta[:, _BETA_IDX]
    se_c = se_all[:, _BETA_IDX]
    z = np.where(se_c > 0, est / se_c, np.nan)
    p = 2.0 * norm.sf(np.abs(z))

    return CommonFactorGWASResult(
        est=est,
        se_c=se_c,
        z=z,
        p=p,
        fail=fail.astype(np.int8),
        warning=warning,
        n_iter=n_iter,
    )
