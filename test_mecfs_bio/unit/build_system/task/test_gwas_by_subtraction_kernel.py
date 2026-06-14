"""
Pure-numpy math tests for the GWAS-by-subtraction kernel. The headline R
comparison via GenomicSEM::userGWAS lives in
test_gwas_by_subtraction_kernel_r_comparison.py.
"""

from __future__ import annotations

import numpy as np
import pytest

from mecfs_bio.build_system.task.r_tasks.genomic_sem._gwas_by_subtraction_kernel import (
    SubtractionLoadings,
    _jacobian_betas_wrt_s,
    _jacobian_betas_wrt_s_finite_difference,
    _solve_betas,
    _solve_betas_from_s,
    _solve_loadings,
    fit_gwas_by_subtraction,
)

# ---- pure-numpy tests -------------------------------------------------------


def test_solve_loadings_identity():
    """When S_LD is diagonal, a_F = 0 and b and a_R match the trait SDs."""
    S_LD = np.diag([0.08, 0.06])
    loadings = _solve_loadings(S_LD)
    assert loadings.b == pytest.approx(np.sqrt(0.06))
    assert loadings.a_F == pytest.approx(0.0)
    assert loadings.a_R == pytest.approx(np.sqrt(0.08))


def test_solve_loadings_with_covariance():
    S_LD = np.array([[0.08, 0.03], [0.03, 0.06]])
    loadings = _solve_loadings(S_LD)
    assert loadings.b == pytest.approx(np.sqrt(0.06))
    assert loadings.a_F == pytest.approx(0.03 / np.sqrt(0.06))
    expected_a_R_sq = 0.08 - (0.03**2 / 0.06)
    assert loadings.a_R == pytest.approx(np.sqrt(expected_a_R_sq))


def test_solve_betas_simple():
    """beta_F = beta_SNP[:,1]/b; beta_R = subtraction formula."""
    loadings = SubtractionLoadings(b=0.3, a_F=0.2, a_R=0.15)
    beta_SNP = np.array([[0.01, 0.006], [-0.02, 0.003]])
    betas = _solve_betas(beta_SNP, loadings)
    np.testing.assert_allclose(betas.beta_F, beta_SNP[:, 1] / loadings.b)
    expected_r = (
        beta_SNP[:, 0] - (loadings.a_F / loadings.b) * beta_SNP[:, 1]
    ) / loadings.a_R
    np.testing.assert_allclose(betas.beta_R, expected_r)


def test_round_trip_from_known_truth():
    """
    Construct S_LD and beta_SNP from known (b, a_F, a_R, beta_F, beta_R),
    then recover the parameters exactly.
    """
    b_true, a_F_true, a_R_true = 0.25, 0.18, 0.20
    beta_F_true = np.array([0.02, -0.01, 0.005])
    beta_R_true = np.array([0.01, 0.03, -0.02])

    S_LD = np.array(
        [
            [a_F_true**2 + a_R_true**2, a_F_true * b_true],
            [a_F_true * b_true, b_true**2],
        ]
    )
    beta_SNP = np.column_stack(
        [
            a_F_true * beta_F_true + a_R_true * beta_R_true,
            b_true * beta_F_true,
        ]
    )

    loadings = _solve_loadings(S_LD)
    assert loadings.b == pytest.approx(b_true)
    assert loadings.a_F == pytest.approx(a_F_true)
    assert loadings.a_R == pytest.approx(a_R_true)

    betas = _solve_betas(beta_SNP, loadings)
    np.testing.assert_allclose(betas.beta_F, beta_F_true, atol=1e-14)
    np.testing.assert_allclose(betas.beta_R, beta_R_true, atol=1e-14)


def test_fit_returns_finite_se():
    """Full fit on clean synthetic data produces finite, positive SEs."""
    rng = np.random.default_rng(0)
    N = 50
    S_LD = np.array([[0.08, 0.03], [0.03, 0.06]])
    V_LD = np.diag([3e-5, 1e-5, 4e-5])
    I_LD = np.eye(2) * 1.05
    varSNP = rng.uniform(0.1, 0.5, N)
    beta_SNP = rng.normal(0, 0.01, (N, 2))
    SE_SNP = np.full((N, 2), 0.01)

    result = fit_gwas_by_subtraction(
        S_LD=S_LD,
        V_LD=V_LD,
        I_LD=I_LD,
        beta_SNP=beta_SNP,
        SE_SNP=SE_SNP,
        varSNP=varSNP,
    )
    assert not result.fail.any()
    assert np.all(np.isfinite(result.se_c_F))
    assert np.all(np.isfinite(result.se_c_R))
    assert np.all(result.se_c_F > 0)
    assert np.all(result.se_c_R > 0)
    assert np.all(np.isfinite(result.n_eff_F))
    assert np.all(result.n_eff_F > 0)


def test_n_eff_smaller_for_remainder_factor():
    """
    The remainder factor (R) has larger SE than the common factor (F) for the
    same SNP, so N_eff(R) < N_eff(F) in general.
    """
    rng = np.random.default_rng(42)
    N = 200
    S_LD = np.array([[0.08, 0.04], [0.04, 0.06]])
    V_LD = np.diag([3e-5, 1e-5, 4e-5])
    I_LD = np.eye(2) * 1.0
    varSNP = rng.uniform(0.1, 0.5, N)
    beta_SNP = rng.normal(0, 0.01, (N, 2))
    SE_SNP = np.full((N, 2), 0.01)

    result = fit_gwas_by_subtraction(
        S_LD=S_LD,
        V_LD=V_LD,
        I_LD=I_LD,
        beta_SNP=beta_SNP,
        SE_SNP=SE_SNP,
        varSNP=varSNP,
    )
    assert np.median(result.n_eff_R) < np.median(result.n_eff_F)


# ---- analytic vs finite-difference Jacobian --------------------------------


def test_analytic_jacobian_matches_finite_difference():
    """Analytic d(beta)/d(s) agrees with centred finite differences at an
    ordinary interior point."""
    rng = np.random.default_rng(7)
    N = 40
    varSNP = rng.uniform(0.1, 0.5, N)
    beta_SNP = rng.normal(0, 0.01, (N, 2))
    s2 = varSNP * beta_SNP[:, 0]
    s3 = varSNP * beta_SNP[:, 1]
    s4, s5, s6 = 0.08, 0.03, 0.06  # a_R^2 = 0.08 - 0.03^2/0.06 = 0.065

    G = _jacobian_betas_wrt_s(s2, s3, s4, s5, s6, varSNP)
    G_fd = _jacobian_betas_wrt_s_finite_difference(s2, s3, s4, s5, s6, varSNP)
    np.testing.assert_allclose(G, G_fd, rtol=1e-5, atol=1e-8)


def test_analytic_jacobian_stays_correct_for_small_a_r():
    """
    When a_R is small (Var(T1) only slightly above Cov^2/Var(T2)), centred
    finite differences over s4/s5/s6 can straddle the `max(..., 1e-30)` clamp
    in `_solve_betas_from_s` and degrade. The analytic Jacobian does not. We
    verify the analytic form against a very small symmetric step that does not
    cross the clamp.
    """
    rng = np.random.default_rng(11)
    N = 25
    varSNP = rng.uniform(0.1, 0.5, N)
    beta_SNP = rng.normal(0, 0.01, (N, 2))
    s2 = varSNP * beta_SNP[:, 0]
    s3 = varSNP * beta_SNP[:, 1]
    # a_R^2 = s4 - s5^2/s6 = 0.0610 - 0.0600 = 1e-3  -> a_R ~ 0.0316 (small).
    s4, s5, s6 = 0.0610, 0.06, 0.06

    G = _jacobian_betas_wrt_s(s2, s3, s4, s5, s6, varSNP)

    def centred_col(idx: int, h: float) -> np.ndarray:
        """Tiny centred difference for the idx-th element of s, staying well
        inside a_R^2 > 0 so the clamp is never hit."""
        p2, p3, p4, p5, p6 = s2, s3, s4, s5, s6
        m2, m3, m4, m5, m6 = s2, s3, s4, s5, s6
        if idx == 0:
            p2, m2 = s2 + h, s2 - h
        elif idx == 1:
            p3, m3 = s3 + h, s3 - h
        elif idx == 2:
            p4, m4 = s4 + h, s4 - h
        elif idx == 3:
            p5, m5 = s5 + h, s5 - h
        else:
            p6, m6 = s6 + h, s6 - h
        plus = _solve_betas_from_s(p2, p3, p4, p5, p6, varSNP)
        minus = _solve_betas_from_s(m2, m3, m4, m5, m6, varSNP)
        return np.stack(
            [
                (plus.beta_F - minus.beta_F) / (2 * h),
                (plus.beta_R - minus.beta_R) / (2 * h),
            ]
        )

    steps = [1e-7, 1e-7, 1e-9, 1e-9, 1e-9]
    for j, h in enumerate(steps):
        ref = centred_col(j, h)
        np.testing.assert_allclose(G[:, :, j].T, ref, rtol=1e-3, atol=1e-6)
