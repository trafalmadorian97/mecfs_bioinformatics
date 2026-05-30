"""
Tests for the numpy GWAS-by-subtraction kernel.

Pure-numpy math tests plus a headline R comparison via GenomicSEM::userGWAS
with the Cholesky model.
"""

from __future__ import annotations

import numpy as np
import pytest

from mecfs_bio.build_system.task.r_tasks.genomic_sem._gwas_by_subtraction_kernel import (
    _jacobian_betas_wrt_s,
    _jacobian_betas_wrt_s_fd,
    _solve_betas,
    _solve_betas_from_s,
    _solve_loadings,
    fit_gwas_by_subtraction,
)

# ---- pure-numpy tests -------------------------------------------------------


def test_solve_loadings_identity():
    """When S_LD is diagonal, a_C = 0 and b and a_NC match the trait SDs."""
    S_LD = np.diag([0.08, 0.06])
    b, a_C, a_NC = _solve_loadings(S_LD)
    assert b == pytest.approx(np.sqrt(0.06))
    assert a_C == pytest.approx(0.0)
    assert a_NC == pytest.approx(np.sqrt(0.08))


def test_solve_loadings_with_covariance():
    S_LD = np.array([[0.08, 0.03], [0.03, 0.06]])
    b, a_C, a_NC = _solve_loadings(S_LD)
    assert b == pytest.approx(np.sqrt(0.06))
    assert a_C == pytest.approx(0.03 / np.sqrt(0.06))
    expected_a_NC_sq = 0.08 - (0.03**2 / 0.06)
    assert a_NC == pytest.approx(np.sqrt(expected_a_NC_sq))


def test_solve_betas_simple():
    """beta_C = beta_SNP[:,1]/b; beta_NC = subtraction formula."""
    b, a_C, a_NC = 0.3, 0.2, 0.15
    beta_SNP = np.array([[0.01, 0.006], [-0.02, 0.003]])
    beta_C, beta_NC = _solve_betas(beta_SNP, b, a_C, a_NC)
    np.testing.assert_allclose(beta_C, beta_SNP[:, 1] / b)
    expected_nc = (beta_SNP[:, 0] - (a_C / b) * beta_SNP[:, 1]) / a_NC
    np.testing.assert_allclose(beta_NC, expected_nc)


def test_round_trip_from_known_truth():
    """
    Construct S_LD and beta_SNP from known (b, a_C, a_NC, beta_C, beta_NC),
    then recover the parameters exactly.
    """
    b_true, a_C_true, a_NC_true = 0.25, 0.18, 0.20
    beta_C_true = np.array([0.02, -0.01, 0.005])
    beta_NC_true = np.array([0.01, 0.03, -0.02])

    S_LD = np.array(
        [
            [a_C_true**2 + a_NC_true**2, a_C_true * b_true],
            [a_C_true * b_true, b_true**2],
        ]
    )
    beta_SNP = np.column_stack(
        [
            a_C_true * beta_C_true + a_NC_true * beta_NC_true,
            b_true * beta_C_true,
        ]
    )

    b, a_C, a_NC = _solve_loadings(S_LD)
    assert b == pytest.approx(b_true)
    assert a_C == pytest.approx(a_C_true)
    assert a_NC == pytest.approx(a_NC_true)

    beta_C, beta_NC = _solve_betas(beta_SNP, b, a_C, a_NC)
    np.testing.assert_allclose(beta_C, beta_C_true, atol=1e-14)
    np.testing.assert_allclose(beta_NC, beta_NC_true, atol=1e-14)


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
    assert np.all(result.fail == 0)
    assert np.all(np.isfinite(result.se_c_C))
    assert np.all(np.isfinite(result.se_c_NC))
    assert np.all(result.se_c_C > 0)
    assert np.all(result.se_c_NC > 0)
    assert np.all(np.isfinite(result.n_eff_C))
    assert np.all(result.n_eff_C > 0)


def test_n_eff_smaller_for_subtracted_factor():
    """
    The subtracted factor (NC) has larger SE than the reference-aligned
    factor (C) for the same SNP, so N_eff(NC) < N_eff(C) in general.
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
    assert np.median(result.n_eff_NC) < np.median(result.n_eff_C)


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
    s4, s5, s6 = 0.08, 0.03, 0.06  # a_NC^2 = 0.08 - 0.03^2/0.06 = 0.065

    G = _jacobian_betas_wrt_s(s2, s3, s4, s5, s6, varSNP)
    G_fd = _jacobian_betas_wrt_s_fd(s2, s3, s4, s5, s6, varSNP)
    np.testing.assert_allclose(G, G_fd, rtol=1e-5, atol=1e-8)


def test_analytic_jacobian_stays_correct_for_small_a_nc():
    """
    When a_NC is small (Var(T1) only slightly above Cov^2/Var(T2)), centred
    finite differences over s4/s5/s6 can straddle the `max(..., 1e-30)` clamp
    in `_solve_betas_from_s` and degrade. The analytic Jacobian does not. We
    verify the analytic form against a tight one-sided complex-free check:
    a very small symmetric step that does not cross the clamp.
    """
    rng = np.random.default_rng(11)
    N = 25
    varSNP = rng.uniform(0.1, 0.5, N)
    beta_SNP = rng.normal(0, 0.01, (N, 2))
    s2 = varSNP * beta_SNP[:, 0]
    s3 = varSNP * beta_SNP[:, 1]
    # a_NC^2 = s4 - s5^2/s6 = 0.0610 - 0.0600 = 1e-3  -> a_NC ~ 0.0316 (small).
    s4, s5, s6 = 0.0610, 0.06, 0.06

    G = _jacobian_betas_wrt_s(s2, s3, s4, s5, s6, varSNP)

    def centred_col(idx: int, h: float) -> np.ndarray:
        """Tiny centred difference for the idx-th element of s, staying well
        inside a_NC^2 > 0 so the clamp is never hit."""
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
        bc_p, bnc_p = _solve_betas_from_s(p2, p3, p4, p5, p6, varSNP)
        bc_m, bnc_m = _solve_betas_from_s(m2, m3, m4, m5, m6, varSNP)
        return np.stack([(bc_p - bc_m) / (2 * h), (bnc_p - bnc_m) / (2 * h)])

    steps = [1e-7, 1e-7, 1e-9, 1e-9, 1e-9]
    for j, h in enumerate(steps):
        ref = centred_col(j, h)
        np.testing.assert_allclose(G[:, :, j].T, ref, rtol=1e-3, atol=1e-6)


# ---- R comparison test -------------------------------------------------------


def _have_rpy2_genomic_sem() -> bool:
    try:
        import rpy2.robjects as ro  # noqa: F401
        from rpy2.robjects.packages import importr

        importr("GenomicSEM")
        return True
    except Exception:
        return False


@pytest.mark.skipif(
    not _have_rpy2_genomic_sem(), reason="rpy2/GenomicSEM not available"
)
def test_matches_user_gwas_cholesky_on_small_slice():
    """
    Run GenomicSEM::userGWAS with the Cholesky subtraction model on
    synthetic data, and assert our kernel matches for both C~SNP and NC~SNP.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    gsem = importr("GenomicSEM")
    base = importr("base")

    rng = np.random.default_rng(42)
    N_SNPS = 30

    # Synthetic LDSC output for 2 traits.
    S_LD = np.array([[0.08, 0.03], [0.03, 0.06]])
    V_LD = np.diag([3e-5, 1e-5, 4e-5])
    I_LD = np.eye(2) * 1.05

    maf = rng.uniform(0.05, 0.5, N_SNPS)
    varSNP = 2.0 * maf * (1.0 - maf)
    beta_SNP = rng.normal(0, 0.01, (N_SNPS, 2))
    SE_SNP = np.full((N_SNPS, 2), 0.01)

    # ---- Build R inputs ----
    trait_names = ["t1", "t2"]
    r_S = base.matrix(ro.FloatVector(S_LD.flatten("F")), nrow=2)
    r_S.rownames = ro.StrVector(trait_names)
    r_S.colnames = ro.StrVector(trait_names)
    r_V = base.matrix(ro.FloatVector(V_LD.flatten("F")), nrow=3)
    r_I = base.matrix(ro.FloatVector(I_LD.flatten("F")), nrow=2)
    r_I.rownames = ro.StrVector(trait_names)
    r_I.colnames = ro.StrVector(trait_names)
    covstruc = ro.ListVector({"V": r_V, "S": r_S, "I": r_I})

    import pandas as pd

    snps_df = pd.DataFrame(
        {
            "SNP": [f"rs{i}" for i in range(N_SNPS)],
            "CHR": rng.integers(1, 23, N_SNPS),
            "BP": rng.integers(1_000_000, 100_000_000, N_SNPS),
            "MAF": maf,
            "A1": ["A"] * N_SNPS,
            "A2": ["G"] * N_SNPS,
            "beta.t1": beta_SNP[:, 0],
            "se.t1": SE_SNP[:, 0],
            "beta.t2": beta_SNP[:, 1],
            "se.t2": SE_SNP[:, 1],
        }
    )
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_snps = ro.conversion.get_conversion().py2rpy(snps_df)

    # Cholesky subtraction model in lavaan syntax.
    model = (
        "C  =~ NA*t1 + t2\n"
        "NC =~ NA*t1\n"
        "C  ~ SNP\n"
        "NC ~ SNP\n"
        "NC ~~ 1*NC\n"
        "C  ~~ 1*C\n"
        "C  ~~ 0*NC\n"
        "t2 ~~ 0*t1\n"
        "t2 ~~ 0*t2\n"
        "t1 ~~ 0*t1\n"
        "SNP ~~ SNP\n"
    )
    r_result = gsem.userGWAS(
        covstruc=covstruc,
        SNPs=r_snps,
        estimation="DWLS",
        model=model,
        sub=ro.StrVector(["C~SNP", "NC~SNP"]),
        parallel=False,
        toler=1e-50,
    )

    # userGWAS returns an R list with one df per sub component.
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_c_df = ro.conversion.get_conversion().rpy2py(r_result.rx2(1))
        r_nc_df = ro.conversion.get_conversion().rpy2py(r_result.rx2(2))

    # userGWAS names the sandwich SE column "SE" (not "se_c" like commonfactorGWAS).
    r_est_C = np.asarray(r_c_df["est"], dtype=float)
    r_se_c_C = np.asarray(r_c_df["SE"], dtype=float)
    r_est_NC = np.asarray(r_nc_df["est"], dtype=float)
    r_se_c_NC = np.asarray(r_nc_df["SE"], dtype=float)

    # ---- Run Python kernel ----
    py = fit_gwas_by_subtraction(
        S_LD=S_LD,
        V_LD=V_LD,
        I_LD=I_LD,
        beta_SNP=beta_SNP,
        SE_SNP=SE_SNP,
        varSNP=varSNP,
    )

    # ---- Compare ----
    ok = (py.fail == 0) & np.isfinite(r_est_C) & np.isfinite(r_est_NC)
    assert ok.sum() >= 0.9 * N_SNPS, f"Too few comparable: {ok.sum()}/{N_SNPS}"

    # userGWAS with fix_measurement=TRUE obtains baseline loadings via lavaan's
    # DWLS optimizer; our kernel uses the equivalent closed-form solution.
    # Both are correct at exact identification, but numerical paths differ
    # slightly, giving ~0.1% discrepancy in the derived betas and SEs.
    np.testing.assert_allclose(py.beta_C[ok], r_est_C[ok], atol=5e-5, rtol=2e-3)
    np.testing.assert_allclose(py.beta_NC[ok], r_est_NC[ok], atol=5e-5, rtol=2e-3)
    np.testing.assert_allclose(py.se_c_C[ok], r_se_c_C[ok], atol=5e-4, rtol=5e-2)
    np.testing.assert_allclose(py.se_c_NC[ok], r_se_c_NC[ok], atol=5e-4, rtol=5e-2)
