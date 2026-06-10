"""
Vectorized GWAS-by-subtraction via Cholesky decomposition.

Given GWAS summary statistics for two traits, this module decomposes each
SNP's genetic effect into a component shared with the reference trait and a
residual component unique to the composite trait.

We name the two latent factors:
- F ("factor"): the genetic factor common to both traits. It is defined by
  the reference trait T2, which is a pure indicator of F.
- R ("remainder"): the residual genetic factor carried only by the composite
  trait T1, orthogonal to F.

In the canonical application (Demange et al. 2021, Nat Genet), T1 =
educational attainment, T2 = cognitive performance, F is the "cognitive"
component and R the "non-cognitive" component. The F/R naming is
use-case-neutral so the same task serves other trait pairs.

References:
- Demange et al. 2021, Nat Genet 53:35-44
- GenomicSEM wiki, section 5 (userGWAS with Cholesky model)
- https://rpubs.com/michelnivard/565885 (tutorial)


# Key matrices and arrays

N = number of SNPs. Traits are ordered (T1 = composite, T2 = reference).

    Symbol      Shape       Meaning
    --------    --------    --------------------------------------------------
    S_LD        (2, 2)      genetic covariance matrix from LDSC (rows/cols =
                            traits; S_LD[0,0]=Var(T1), S_LD[1,1]=Var(T2))
    V_LD        (3, 3)      sampling covariance of vech(S_LD) from LDSC, over
                            (S_LD[0,0], S_LD[0,1], S_LD[1,1])
    I_LD        (2, 2)      LDSC intercept matrix (cross-trait sample overlap)
    V_SNP       (N, 2, 2)   per-SNP sampling covariance of
                            (Cov(SNP,T1), Cov(SNP,T2))
    V_s         (N, 5, 5)   blockdiag(V_SNP, V_LD): full sampling covariance
                            of the moment vector s
    beta_SNP    (N, 2)      per-SNP marginal GWAS betas (col 0 = T1, col 1 = T2)
    SE_SNP      (N, 2)      per-SNP standard errors
    varSNP      (N,)        2 * MAF * (1 - MAF)


# The structural equation model

Two latent factors F (common) and R (remainder), both standardised to unit
variance and constrained to be uncorrelated:

    T2 = b    * F                               (T2 is a pure indicator of F)
    T1 = a_F  * F  +  a_R * R                   (T1 loads on both factors)

    Var(F) = 1,  Var(R) = 1,  Cov(F, R) = 0
    No residual variance in T1 or T2:  all genetic variance is explained
    by the two factors.

Adding a SNP regressor:

    F = beta_F * SNP + zeta_F                    Var(zeta_F) = 1 - beta_F^2 * varSNP
    R = beta_R * SNP + zeta_R                    Var(zeta_R) = 1 - beta_R^2 * varSNP


# Closed-form solution

Because the model is exactly identified (5 informative moment equations,
5 free parameters), the estimates are explicit:

    b       = sqrt(S_LD[1,1])
    a_F     = S_LD[0,1] / b
    a_R     = sqrt(S_LD[0,0] - a_F^2)

    beta_F  = beta_SNP[:, 1] / b                          (reference-aligned)
    beta_R  = (beta_SNP[:, 0] - (a_F/b) * beta_SNP[:, 1]) / a_R   (subtraction)

No iterative optimisation is needed — the entire genome can be processed in
one vectorised pass.


# Standard errors via the delta method

The parameters theta = (b, a_F, a_R, beta_F, beta_R) are smooth functions
of the data s = (s2, s3, s4, s5, s6) = (Cov(SNP,T1), Cov(SNP,T2),
Var(T1), Cov(T1,T2), Var(T2)). The delta method gives:

    Var(theta_hat) = G * V_s * G'

where G = d theta / d s is the (5 x 5) Jacobian evaluated at the data, and
V_s = blockdiag(V_SNP, V_LD) is the sampling covariance of s. We report SE
for beta_F and beta_R from the diagonal of the relevant 2x2 submatrix.

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
from attrs import frozen
from scipy.stats import norm


@frozen
class SubtractionLoadings:
    """
    Cholesky factor loadings (constant across SNPs).

    b   : loading of the reference trait (T2) on the common factor F.
    a_F : loading of the composite trait (T1) on the common factor F.
    a_R : loading of the composite trait (T1) on the remainder factor R.


    Derived by equating the theoretical and empirical 2x2 covariance matricies

    [a_F^2 + a_R^2  a_F*b]        ==      [S_{LD,1,1}        S_{LD,1,2}]
    [a_F*b          b^2  ]                [S_{LD,1,2}        S_{LD,2,2}]

    Which yields
    b=\sqrt{S_{LD,2,2}}
    a_F=S_{LD,1,2}/\sqrt{S_{LD,2,2}}
    a_R=\sqrt{S_{LD,1,1}-\frac{S_{LD,1,2}^2}{S_{LD,2,2}}}

    """

    b: float
    a_F: float
    a_R: float

    def __attrs_post_init__(self):
        assert self.b>=0, f"b is equal to the variance of trait 2, and so must be positive.  Got b={self.b}"
        assert self.a_R>=0, f"a_R must be positive.  If it is negative, this implies the genetic covariance matrix is not positive definite."


@frozen
class FactorBetas:
    """Per-SNP factor effects (each shape (N,))."""

    beta_F: np.ndarray  # effect on the common factor F
    beta_R: np.ndarray  # effect on the remainder factor R

    def __attrs_post_init__(self):
        assert self.beta_F.ndim==1
        assert (
            self.beta_F.shape == self.beta_R.shape
        )


@frozen
class GWASBySubtractionResult:
    """Per-SNP results from `fit_gwas_by_subtraction`."""

    loadings: SubtractionLoadings  # the Cholesky loadings used for every SNP

    # Per-SNP arrays (length N).
    beta_F: np.ndarray  # common-factor effect
    se_c_F: np.ndarray  # sandwich SE for beta_F
    z_F: np.ndarray # Z scores of the beta-F
    p_F: np.ndarray  # p-values of the beta-Fs
    beta_R: np.ndarray  # remainder-factor effect
    se_c_R: np.ndarray  # sandwich SE for beta_R
    z_R: np.ndarray  # z scores of the beta_Rs
    p_R: np.ndarray  # p values of the beta_Rs
    n_eff_F: np.ndarray  # effective sample size for beta_F
    n_eff_R: np.ndarray  # effective sample size for beta_R
    fail: np.ndarray  # bool, True if SE was non-finite/non-positive for this SNP

    def __attrs_post_init__(self):
        assert (
            self.beta_F.ndim==1
        )
        assert (self.beta_F.shape
                ==self.se_c_F.shape
                ==self.z_F.shape
                ==self.p_F.shape
                ==self.beta_R.shape
                ==self.se_c_R.shape
                ==self.z_R.shape
                ==self.p_R.shape
                ==self.n_eff_F.shape
                ==self.n_eff_R.shape
                ==self.fail.shape)


def _solve_loadings(S_LD: np.ndarray) -> SubtractionLoadings:
    """
    Compute the Cholesky loadings from the 2x2 genetic covariance matrix.

    Convention: source[0] = composite (T1), source[1] = reference (T2).
    Factor F is shared between both traits; factor R is the T1 residual.



    Derived by equating the theoretical and empirical 2x2 covariance matricies

    [a_F^2 + a_R^2  a_F*b]        ==      [S_{LD,1,1}        S_{LD,1,2}]
    [a_F*b          b^2  ]                [S_{LD,1,2}        S_{LD,2,2}]

    Which yields
    b=\sqrt{S_{LD,2,2}}
    a_F=S_{LD,1,2}/\sqrt{S_{LD,2,2}}
    a_R=\sqrt{S_{LD,1,1}-\frac{S_{LD,1,2}^2}{S_{LD,2,2}}}
    """
    assert S_LD.shape == (2, 2), f"S_LD must be (2, 2); got {S_LD.shape}"
    s44 = float(S_LD[0, 0])  # Var(T1)
    s45 = float(S_LD[0, 1])  # Cov(T1, T2)
    s55 = float(S_LD[1, 1])  # Var(T2)

    b = np.sqrt(s55)
    a_F = s45 / b
    a_R_sq = s44 - a_F**2
    assert a_R_sq > 0, (
        f"S_LD implies a_R^2 = {a_R_sq:.6e} <= 0; the Cholesky "
        f"decomposition requires Var(T1) > Cov(T1,T2)^2 / Var(T2)."
    )
    a_R = np.sqrt(a_R_sq)
    return SubtractionLoadings(b=b, a_F=a_F, a_R=a_R)


def _solve_betas(beta_SNP: np.ndarray, loadings: SubtractionLoadings) -> FactorBetas:
    """
    Per-SNP closed-form factor effects.

    beta_SNP : (N, 2) with column 0 = T1 marginal beta, column 1 = T2 marginal beta.
    Returns FactorBetas with beta_F, beta_R each shape (N,).


    The formulae used here are derived by equating the first columns of the empirical and theoretical
    3x3 covariance matrices

    That is:

    [
    Var_{SNP}
    Var_{SNP} * (a_F * \beta_{F,i} + a_R * \beta_{R,i} )
    Var_{SNP} * b * \beta_{F,i}
    ]

    ==

    [
    Var_{SNP}
    \beta_{T_1, i} * Var_{SNP}
    \beta_{T_2, i} * Var_{SNP}
    ].


    Solving these yields

    \beta_{F,i}  = \beta_{ T_2,i } / b
    \beta_{R,i} = 1/a_R * (\beta_{ T_1,i } -a_F * \beta_{T_2,i}  / b )


    """
    assert beta_SNP.ndim == 2 and beta_SNP.shape[1] == 2, (
        f"beta_SNP must be (N, 2); got {beta_SNP.shape}"
    )
    b, a_F, a_R = loadings.b, loadings.a_F, loadings.a_R
    beta_F = beta_SNP[:, 1] / b
    beta_R = (beta_SNP[:, 0] - (a_F / b) * beta_SNP[:, 1]) / a_R
    return FactorBetas(beta_F=beta_F, beta_R=beta_R)


def _solve_betas_from_s(
    s2: np.ndarray,
    s3: np.ndarray,
    s4: float,
    s5: float,
    s6: float,
    varSNP: np.ndarray,
) -> FactorBetas:
    """
    Solve for (beta_F, beta_R) from the raw moment vector s.

    s2 = Cov(SNP, T1) = varSNP * beta_SNP[:,0]   (N,)
    s3 = Cov(SNP, T2) = varSNP * beta_SNP[:,1]   (N,)
    s4 = S_LD[0,0],  s5 = S_LD[0,1],  s6 = S_LD[1,1]   (scalars)

    Differs from _solve_betas, in that we solve for a and b jointly with solving for the betas

    Used only in finite difference Jacobian computation
    """
    assert s2.shape == s3.shape == varSNP.shape, (
        f"s2, s3, varSNP must share shape; got {s2.shape}, {s3.shape}, {varSNP.shape}"
    )
    b = np.sqrt(s6)
    a_F = s5 / b
    a_R = np.sqrt(max(s4 - a_F**2, 1e-30))
    beta_F = s3 / (b * varSNP)
    beta_R = (s2 / varSNP - (a_F / b) * s3 / varSNP) / a_R
    return FactorBetas(beta_F=beta_F, beta_R=beta_R)


def _jacobian_betas_wrt_s(
    s2: np.ndarray,
    s3: np.ndarray,
    s4: float,
    s5: float,
    s6: float,
    varSNP: np.ndarray,
) -> np.ndarray:
    """
    Analytic d(beta_F, beta_R)/d(s2, s3, s4, s5, s6). Returns shape (N, 2, 5).

    The estimates are elementary functions of the data:

        b      = sqrt(s6)
        a_F    = s5 / b               (so a_F / b = s5 / s6)
        a_R    = sqrt(s4 - s5^2 / s6)
        beta_F  = s3 / (b * v)
        beta_R  = (s2 / v - (s5 / s6) * s3 / v) / a_R

    with v = varSNP. Differentiating directly (v is known data, not part of
    s) gives the partials below; see the module test for a finite-difference
    cross-check. Unlike finite differences, this never evaluates the
    `max(s4 - a_F^2, ...)` clamp inside `_solve_betas_from_s`, so it stays
    correct even when a_R is small.
    """
    assert s2.shape == s3.shape == varSNP.shape, (
        f"s2, s3, varSNP must share shape; got {s2.shape}, {s3.shape}, {varSNP.shape}"
    )
    N = s2.shape[0]
    v = varSNP
    b = np.sqrt(s6)
    a_R_sq = s4 - s5**2 / s6
    a_R = np.sqrt(a_R_sq)
    beta_R = (s2 / v - (s5 / s6) * s3 / v) / a_R

    G = np.zeros((N, 2, 5))

    # beta_F = s3 / (sqrt(s6) * v): depends only on s3 and s6.
    beta_F = s3 / (b * v)
    G[:, 0, 1] = 1.0 / (b * v)  # d/ds3
    G[:, 0, 4] = -beta_F / (2.0 * s6)  # d/ds6

    # beta_R partials.
    G[:, 1, 0] = 1.0 / (v * a_R)  # d/ds2
    G[:, 1, 1] = -s5 / (s6 * v * a_R)  # d/ds3
    G[:, 1, 2] = -beta_R / (2.0 * a_R_sq)  # d/ds4
    G[:, 1, 3] = -s3 / (s6 * v * a_R) + beta_R * s5 / (s6 * a_R_sq)  # d/ds5
    G[:, 1, 4] = s5 * s3 / (s6**2 * v * a_R) - beta_R * s5**2 / (
        2.0 * s6**2 * a_R_sq
    )  # d/ds6
    return G


def _jacobian_betas_wrt_s_finite_difference(
    s2: np.ndarray,
    s3: np.ndarray,
    s4: float,
    s5: float,
    s6: float,
    varSNP: np.ndarray,
) -> np.ndarray:
    """
    Finite-difference reference for `_jacobian_betas_wrt_s`, kept as a test
    oracle only (not used in production). Centred differences with a
    scale-aware step. Returns shape (N, 2, 5).
    """
    assert s2.shape == s3.shape == varSNP.shape, (
        f"s2, s3, varSNP must share shape; got {s2.shape}, {s3.shape}, {varSNP.shape}"
    )
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
        plus = _solve_betas_from_s(d2, d3, d4, d5, d6, varSNP)
        m2 = s2 - h if idx == 0 else s2
        m3 = s3 - h if idx == 1 else s3
        m4 = s4 - h if idx == 2 else s4
        m5 = s5 - h if idx == 3 else s5
        m6 = s6 - h if idx == 4 else s6
        minus = _solve_betas_from_s(m2, m3, m4, m5, m6, varSNP)
        return plus.beta_F, plus.beta_R, minus.beta_F, minus.beta_R

    scales = [
        float(np.median(np.abs(s2))),
        float(np.median(np.abs(s3))),
        abs(s4),
        abs(s5),
        abs(s6),
    ]
    for j in range(5):
        h = max(scales[j] * 1e-5, 1e-8)
        bf_p, br_p, bf_m, br_m = _perturb(j, h)
        G[:, 0, j] = (bf_p - bf_m) / (2 * h)
        G[:, 1, j] = (br_p - br_m) / (2 * h)
    return G


def compute_v_snp_batch(
    SE_SNP: np.ndarray, I_LD: np.ndarray, varSNP: np.ndarray
) -> np.ndarray:
    """
    Per-SNP sampling covariance of (Cov(SNP,T1), Cov(SNP,T2)) for
    GC="standard". Same formula as in `_common_factor_kernel`.

    SE_SNP : (N, 2),  I_LD : (2, 2),  varSNP : (N,).  Returns (N, 2, 2).

    Inputs:

    SE_SNP: (N,2) standard errors of the two GWAS traits T1 and T2
    I_LD: (2,2) LDSC intercept matrix
    varSNP: (N,) Per SNP variances




    """
    assert SE_SNP.ndim == 2 and SE_SNP.shape[1] == 2, (
        f"SE_SNP must be (N, 2); got {SE_SNP.shape}"
    )
    assert I_LD.shape == (2, 2), f"I_LD must be (2, 2); got {I_LD.shape}"
    assert varSNP.shape == (SE_SNP.shape[0],), (
        f"varSNP must be (N,); got {varSNP.shape}"
    )
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
) -> GWASBySubtractionResult:
    """
    GWAS-by-subtraction for k=2 traits via Cholesky decomposition.

    Convention: source 0 = composite trait (T1), source 1 = reference trait (T2).
    Factor F is shared between both traits; factor R is the T1 residual
    orthogonal to F. See the module docstring for the matrix shapes.

    Inputs:
        S_LD       : (2, 2) genetic covariance from LDSC
        V_LD       : (3, 3) sampling covariance of vech(S_LD) from LDSC
        I_LD       : (2, 2) LDSC intercepts
        beta_SNP   : (N, 2) per-SNP marginal betas (col 0 = T1, col 1 = T2)
        SE_SNP     : (N, 2) per-SNP standard errors
        varSNP     : (N,) = 2 * MAF * (1 - MAF)

    Returns GWASBySubtractionResult.

    Raises ValueError if the sampling covariance V_s is not positive-definite
    (a malformed-input condition; see below).
    """
    assert beta_SNP.ndim == 2, f"beta_SNP must be 2-D (N, 2); got {beta_SNP.shape}"
    N, k = beta_SNP.shape
    assert k == 2, f"GWAS-by-subtraction requires exactly 2 traits; got {k}"
    assert S_LD.shape == (2, 2), f"S_LD must be (2, 2); got {S_LD.shape}"
    assert V_LD.shape == (3, 3), f"V_LD must be (3, 3); got {V_LD.shape}"
    assert I_LD.shape == (2, 2), f"I_LD must be (2, 2); got {I_LD.shape}"
    assert SE_SNP.shape == (N, 2), f"SE_SNP must be {(N, 2)}; got {SE_SNP.shape}"
    assert varSNP.shape == (N,), f"varSNP must be {(N,)}; got {varSNP.shape}"

    # Step 1: loadings from S_LD (constant across SNPs).
    loadings = _solve_loadings(S_LD)

    # Step 2: per-SNP betas (vectorised, no loop).
    betas = _solve_betas(beta_SNP, loadings)
    beta_F, beta_R = betas.beta_F, betas.beta_R

    # Step 3: build per-SNP sampling covariance V_s = blockdiag(V_SNP, V_LD).
    V_SNP_batch = compute_v_snp_batch(SE_SNP, I_LD, varSNP)  # (N, 2, 2)
    V_s = np.zeros((N, 5, 5))
    V_s[:, :2, :2] = V_SNP_batch
    V_s[:, 2:, 2:] = V_LD[None, :, :]

    # V_s positive-definiteness is determined by the shared LDSC inputs
    # (V_LD, and whether |I_LD off-diagonal| < 1), not by per-SNP data, so a
    # non-PD V_s is a malformed-input condition affecting every SNP. Fail loudly
    # rather than emitting a per-SNP flag that callers might ignore.
    eig_min = np.linalg.eigvalsh(V_s).min(axis=-1)
    if np.any(eig_min <= 0):
        n_bad = int(np.sum(eig_min <= 0))
        raise ValueError(
            f"Sampling covariance V_s is not positive-definite for {n_bad}/{N} "
            f"SNP(s) (min eigenvalue {float(eig_min.min()):.3e}). This is set by "
            f"the shared LDSC inputs (non-PD V_LD, or |I_LD off-diagonal| >= 1), "
            f"so it indicates malformed inputs rather than a per-SNP condition."
        )

    # Step 4: Jacobian of (beta_F, beta_R) w.r.t. s = (s2, s3, s4, s5, s6).
    s2 = varSNP * beta_SNP[:, 0]
    s3 = varSNP * beta_SNP[:, 1]
    s4 = float(S_LD[0, 0])
    s5 = float(S_LD[0, 1])
    s6 = float(S_LD[1, 1])
    G = _jacobian_betas_wrt_s(s2, s3, s4, s5, s6, varSNP)  # (N, 2, 5)

    # Step 5: delta-method variance  Var(beta) = G V_s G'.
    Var_beta = G @ V_s @ np.swapaxes(G, -2, -1)  # (N, 2, 2)
    var_F = Var_beta[:, 0, 0]
    var_R = Var_beta[:, 1, 1]

    # With V_s positive-definite, var_F/var_R are quadratic forms x'V_s x and
    # are strictly positive for finite, nonzero G. `fail` therefore only flags
    # degenerate per-SNP inputs (non-finite SE, or a monomorphic SNP with
    # varSNP == 0 making G blow up) — it should be empty on clean data.
    fail = ~(np.isfinite(var_F) & np.isfinite(var_R) & (var_F > 0) & (var_R > 0))
    good = ~fail

    se_c_F = np.full(N, np.nan)
    se_c_R = np.full(N, np.nan)
    se_c_F[good] = np.sqrt(var_F[good])
    se_c_R[good] = np.sqrt(var_R[good])

    z_F = beta_F / se_c_F
    z_R = beta_R / se_c_R
    p_F = 2.0 * norm.sf(np.abs(z_F))
    p_R = 2.0 * norm.sf(np.abs(z_R))

    # Effective sample size: N_eff = 1 / (varSNP * se^2).
    n_eff_F = np.full(N, np.nan)
    n_eff_R = np.full(N, np.nan)
    n_eff_F[good] = 1.0 / (varSNP[good] * se_c_F[good] ** 2)
    n_eff_R[good] = 1.0 / (varSNP[good] * se_c_R[good] ** 2)

    return GWASBySubtractionResult(
        loadings=loadings,
        beta_F=beta_F,
        se_c_F=se_c_F,
        z_F=z_F,
        p_F=p_F,
        beta_R=beta_R,
        se_c_R=se_c_R,
        z_R=z_R,
        p_R=p_R,
        n_eff_F=n_eff_F,
        n_eff_R=n_eff_R,
        fail=fail,
    )
