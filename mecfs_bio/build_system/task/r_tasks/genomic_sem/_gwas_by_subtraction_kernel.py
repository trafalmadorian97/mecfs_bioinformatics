"""
Vectorized GWAS-by-subtraction via Cholesky decomposition.

Given GWAS summary statistics for two traits (a "composite" T1 and a
"reference" T2), this module decomposes each SNP's genetic effect into a
component aligned with T2 and a residual component orthogonal to T2. In the
canonical application (Demange et al. 2021, Nat Genet), T1 = educational
attainment, T2 = cognitive performance, and the residual component is
interpreted as "non-cognitive" genetic influence on education.

References:
- Demange et al. 2021, Nat Genet 53:35-44
- GenomicSEM wiki, section 5 (userGWAS with Cholesky model)
- https://rpubs.com/michelnivard/565885 (tutorial)


# The structural equation model

Two latent factors C (reference-aligned) and NC (residual), both
standardised to unit variance and constrained to be uncorrelated:

    T2 = b    * C                               (T2 is a pure indicator of C)
    T1 = a_C  * C  +  a_NC * NC                 (T1 loads on both factors)

    Var(C) = 1,  Var(NC) = 1,  Cov(C, NC) = 0
    No residual variance in T1 or T2:  all genetic variance is explained
    by the two factors.

Adding a SNP regressor:

    C  = beta_C  * SNP + zeta_C                  Var(zeta_C)  = 1 - beta_C^2 * varSNP
    NC = beta_NC * SNP + zeta_NC                  Var(zeta_NC) = 1 - beta_NC^2 * varSNP

The observed data per SNP consists of:
- varSNP = 2 * MAF * (1 - MAF)
- beta_SNP[:, 0], beta_SNP[:, 1]  (marginal GWAS betas for T1 and T2)
- S_LD (2x2 genetic covariance from LDSC)


# Closed-form solution

Because the model is exactly identified (5 informative moment equations,
5 free parameters), the estimates are explicit:

    b       = sqrt(S_LD[1,1])
    a_C     = S_LD[0,1] / b
    a_NC    = sqrt(S_LD[0,0] - a_C^2)

    beta_C  = beta_SNP[:, 1] / b                          (reference-aligned)
    beta_NC = (beta_SNP[:, 0] - (a_C/b) * beta_SNP[:, 1]) / a_NC   (subtraction)

No iterative optimisation is needed — the entire genome can be processed in
one vectorised pass.


# Standard errors via the delta method

The parameters theta = (b, a_C, a_NC, beta_C, beta_NC) are smooth functions
of the data s = (s2, s3, s4, s5, s6) = (Cov(SNP,T1), Cov(SNP,T2),
Var(T1), Cov(T1,T2), Var(T2)). The delta method gives:

    Var(theta_hat) = G * V_s * G'

where G = d theta / d s is the (5 x 5) Jacobian evaluated at the data, and
V_s = blockdiag(V_SNP, V_LD) is the sampling covariance of s (V_SNP per-SNP,
V_LD from LDSC). We report SE for beta_C and beta_NC from the diagonal
of the relevant 2x2 submatrix.

For an exactly-identified model this equals the sandwich SE that
GenomicSEM's userGWAS would produce.


# Effective sample size

For downstream tools (MAGMA, LD score regression on the derived sumstats):

    N_eff_i = 1 / (varSNP_i * se_c_i^2)

This backs out the GWAS-equivalent N from the sandwich SE of the derived
statistic. It will generally be smaller than either input GWAS's N because
the subtraction amplifies noise.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from scipy.stats import norm


class GWASBySubtractionResult(NamedTuple):
    """Per-SNP results from `fit_gwas_by_subtraction`."""

    # Loadings (constant across SNPs, but included for inspection).
    b: float
    a_C: float
    a_NC: float

    # Per-SNP arrays (length N).
    beta_C: np.ndarray  # reference-aligned factor effect
    se_c_C: np.ndarray  # sandwich SE for beta_C
    z_C: np.ndarray
    p_C: np.ndarray
    beta_NC: np.ndarray  # subtracted (residual) factor effect
    se_c_NC: np.ndarray  # sandwich SE for beta_NC
    z_NC: np.ndarray
    p_NC: np.ndarray
    n_eff_C: np.ndarray  # effective sample size for beta_C
    n_eff_NC: np.ndarray  # effective sample size for beta_NC
    fail: np.ndarray  # 1 if SE computation failed for this SNP
    warning: np.ndarray  # 1 if V_Full was non-PD


def _solve_loadings(
    S_LD: np.ndarray,
) -> tuple[float, float, float]:
    """
    Compute the Cholesky loadings from the 2x2 genetic covariance matrix.

    Convention: source[0] = composite (T1), source[1] = reference (T2).
    Factor C aligns with T2; factor NC captures the residual in T1.
    """
    s44 = float(S_LD[0, 0])  # Var(T1)
    s45 = float(S_LD[0, 1])  # Cov(T1, T2)
    s55 = float(S_LD[1, 1])  # Var(T2)

    b = np.sqrt(s55)
    a_C = s45 / b
    a_NC_sq = s44 - a_C**2
    assert a_NC_sq > 0, (
        f"S_LD implies a_NC^2 = {a_NC_sq:.6e} <= 0; the Cholesky "
        f"decomposition requires Var(T1) > Cov(T1,T2)^2 / Var(T2)."
    )
    a_NC = np.sqrt(a_NC_sq)
    return b, a_C, a_NC


def _solve_betas(
    beta_SNP: np.ndarray, b: float, a_C: float, a_NC: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Per-SNP closed-form estimates.

    beta_SNP : (N, 2) with column 0 = T1 marginal beta, column 1 = T2 marginal beta.
    Returns (beta_C, beta_NC), each shape (N,).
    """
    beta_C = beta_SNP[:, 1] / b
    beta_NC = (beta_SNP[:, 0] - (a_C / b) * beta_SNP[:, 1]) / a_NC
    return beta_C, beta_NC


def _solve_betas_from_s(
    s2: np.ndarray,
    s3: np.ndarray,
    s4: float,
    s5: float,
    s6: float,
    varSNP: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Solve for (beta_C, beta_NC) from the raw moment vector s.

    s2 = Cov(SNP, T1) = varSNP * beta_SNP[:,0]
    s3 = Cov(SNP, T2) = varSNP * beta_SNP[:,1]
    s4 = S_LD[0,0],  s5 = S_LD[0,1],  s6 = S_LD[1,1]
    """
    b = np.sqrt(s6)
    a_C = s5 / b
    a_NC = np.sqrt(max(s4 - a_C**2, 1e-30))
    beta_C = s3 / (b * varSNP)
    beta_NC = (s2 / varSNP - (a_C / b) * s3 / varSNP) / a_NC
    return beta_C, beta_NC


def _jacobian_betas_wrt_s(
    s2: np.ndarray,
    s3: np.ndarray,
    s4: float,
    s5: float,
    s6: float,
    varSNP: np.ndarray,
) -> np.ndarray:
    """
    Compute d(beta_C, beta_NC)/d(s2, s3, s4, s5, s6) via centred finite
    differences. Returns shape (N, 2, 5).
    """
    N = s2.shape[0]
    G = np.zeros((N, 2, 5))

    def _perturb(
        idx: int, h: float
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        d2 = s2 + h if idx == 0 else s2
        d3 = s3 + h if idx == 1 else s3
        d4 = s4 + h if idx == 2 else s4
        d5 = s5 + h if idx == 3 else s5
        d6 = s6 + h if idx == 4 else s6
        bc, bnc = _solve_betas_from_s(d2, d3, d4, d5, d6, varSNP)
        m2 = s2 - h if idx == 0 else s2
        m3 = s3 - h if idx == 1 else s3
        m4 = s4 - h if idx == 2 else s4
        m5 = s5 - h if idx == 3 else s5
        m6 = s6 - h if idx == 4 else s6
        bc_m, bnc_m = _solve_betas_from_s(m2, m3, m4, m5, m6, varSNP)
        return bc, bnc, bc_m, bnc_m

    scales = [
        float(np.median(np.abs(s2))),
        float(np.median(np.abs(s3))),
        abs(s4),
        abs(s5),
        abs(s6),
    ]
    for j in range(5):
        h = max(scales[j] * 1e-5, 1e-8)
        bc_p, bnc_p, bc_m, bnc_m = _perturb(j, h)
        G[:, 0, j] = (bc_p - bc_m) / (2 * h)
        G[:, 1, j] = (bnc_p - bnc_m) / (2 * h)
    return G


def compute_v_snp_batch(
    SE_SNP: np.ndarray, I_LD: np.ndarray, varSNP: np.ndarray
) -> np.ndarray:
    """
    Per-SNP sampling covariance of (Cov(SNP,T1), Cov(SNP,T2)) for
    GC="standard". Same formula as in `_common_factor_kernel`.

    Returns (N, 2, 2).
    """
    k = 2
    assert SE_SNP.shape[1] == k
    sqrt_diag_I = np.sqrt(np.diag(I_LD))
    se_scaled = SE_SNP * sqrt_diag_I
    se_outer = se_scaled[:, :, None] * se_scaled[:, None, :]
    intercept_factor = I_LD.copy()
    np.fill_diagonal(intercept_factor, 1.0)
    var_snp_sq = varSNP**2
    return se_outer * intercept_factor[None, :, :] * var_snp_sq[:, None, None]


def fit_gwas_by_subtraction(
    *,
    S_LD: np.ndarray,
    V_LD: np.ndarray,
    I_LD: np.ndarray,
    beta_SNP: np.ndarray,
    SE_SNP: np.ndarray,
    varSNP: np.ndarray,
    varSNPSE2: float = (5e-4) ** 2,
) -> GWASBySubtractionResult:
    """
    GWAS-by-subtraction for k=2 traits via Cholesky decomposition.

    Convention: source 0 = composite trait (T1), source 1 = reference trait (T2).
    Factor C aligns with T2; factor NC captures the T1 residual orthogonal to C.

    Inputs:
        S_LD       : (2, 2) genetic covariance from LDSC
        V_LD       : (3, 3) sampling covariance of vech(S_LD) from LDSC
        I_LD       : (2, 2) LDSC intercepts
        beta_SNP   : (N, 2) per-SNP marginal betas (col 0 = T1, col 1 = T2)
        SE_SNP     : (N, 2) per-SNP standard errors
        varSNP     : (N,) = 2 * MAF * (1 - MAF)
        varSNPSE2  : scalar, variance of the SNP-SNP entry (default (5e-4)^2)

    Returns GWASBySubtractionResult.
    """
    N, k = beta_SNP.shape
    assert k == 2, f"GWAS-by-subtraction requires exactly 2 traits; got {k}"
    assert S_LD.shape == (2, 2)
    assert V_LD.shape == (3, 3)
    assert I_LD.shape == (2, 2)

    # Step 1: loadings from S_LD (constant across SNPs).
    b, a_C, a_NC = _solve_loadings(S_LD)

    # Step 2: per-SNP betas (vectorised, no loop).
    beta_C, beta_NC = _solve_betas(beta_SNP, b, a_C, a_NC)

    # Step 3: build per-SNP sampling covariance V_s = blockdiag(V_SNP, V_LD).
    V_SNP_batch = compute_v_snp_batch(SE_SNP, I_LD, varSNP)  # (N, 2, 2)
    V_s = np.zeros((N, 5, 5))
    V_s[:, :2, :2] = V_SNP_batch
    V_s[:, 2:, 2:] = V_LD[None, :, :]

    # PD check on V_s (diagnostic only).
    eig_min = np.linalg.eigvalsh(V_s).min(axis=-1)
    warning = (eig_min <= 0).astype(np.int8)

    # Step 4: Jacobian of (beta_C, beta_NC) w.r.t. s = (s2, s3, s4, s5, s6).
    s2 = varSNP * beta_SNP[:, 0]
    s3 = varSNP * beta_SNP[:, 1]
    s4 = float(S_LD[0, 0])
    s5 = float(S_LD[0, 1])
    s6 = float(S_LD[1, 1])
    G = _jacobian_betas_wrt_s(s2, s3, s4, s5, s6, varSNP)  # (N, 2, 5)

    # Step 5: delta-method variance  Var(beta) = G V_s G'.
    Var_beta = G @ V_s @ np.swapaxes(G, -2, -1)  # (N, 2, 2)
    var_C = Var_beta[:, 0, 0]
    var_NC = Var_beta[:, 1, 1]

    fail = (
        ~np.isfinite(var_C) | ~np.isfinite(var_NC) | (var_C < 0) | (var_NC < 0)
    ).astype(np.int8)
    se_c_C = np.sqrt(np.where(var_C > 0, var_C, np.nan))
    se_c_NC = np.sqrt(np.where(var_NC > 0, var_NC, np.nan))

    z_C = np.where(se_c_C > 0, beta_C / se_c_C, np.nan)
    z_NC = np.where(se_c_NC > 0, beta_NC / se_c_NC, np.nan)
    p_C = 2.0 * norm.sf(np.abs(z_C))
    p_NC = 2.0 * norm.sf(np.abs(z_NC))

    # Effective sample size: N_eff = 1 / (varSNP * se^2).
    n_eff_C = np.where(se_c_C > 0, 1.0 / (varSNP * se_c_C**2), np.nan)
    n_eff_NC = np.where(se_c_NC > 0, 1.0 / (varSNP * se_c_NC**2), np.nan)

    return GWASBySubtractionResult(
        b=b,
        a_C=a_C,
        a_NC=a_NC,
        beta_C=beta_C,
        se_c_C=se_c_C,
        z_C=z_C,
        p_C=p_C,
        beta_NC=beta_NC,
        se_c_NC=se_c_NC,
        z_NC=z_NC,
        p_NC=p_NC,
        n_eff_C=n_eff_C,
        n_eff_NC=n_eff_NC,
        fail=fail,
        warning=warning,
    )
