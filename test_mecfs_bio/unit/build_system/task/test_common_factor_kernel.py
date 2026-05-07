"""
Tests for the numpy/scipy reimplementation of GenomicSEM::commonfactorGWAS.

Most tests are pure-numpy math sanity checks. The headline test
`test_matches_lavaan_on_small_slice` runs R's commonfactorGWAS via rpy2 on
synthetic data and compares against our kernel.
"""

from __future__ import annotations

import numpy as np
import pytest

from mecfs_bio.build_system.task.r_tasks.genomic_sem._common_factor_kernel import (
    compute_s_fullrun_batch,
    compute_v_full_batch,
    compute_v_snp_batch,
    fit_common_factor_gwas,
    model_implied_vech,
    model_jacobian,
    vech_lower_col_major,
)

# ---- pure-numpy unit tests -------------------------------------------------


def test_vech_lower_col_major_3x3():
    """vech of a 3x3 takes the 6 lower-tri entries in column-major order."""
    M = np.array(
        [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
        ]
    )
    out = vech_lower_col_major(M)
    # Column 0: rows 0,1,2 -> 1, 4, 7
    # Column 1: rows 1,2   -> 5, 8
    # Column 2: rows 2     -> 9
    assert out.tolist() == [1.0, 4.0, 7.0, 5.0, 8.0, 9.0]


def test_vech_lower_col_major_batched():
    M = np.arange(2 * 3 * 3).reshape(2, 3, 3).astype(float)
    out = vech_lower_col_major(M)
    assert out.shape == (2, 6)
    # First batch element is [0..8] reshaped (3, 3); column-major lower tri
    # picks indices (0,0)=0, (1,0)=3, (2,0)=6, (1,1)=4, (2,1)=7, (2,2)=8.
    np.testing.assert_array_equal(out[0], [0, 3, 6, 4, 7, 8])


def test_v_snp_construction_diagonal_only():
    """If I_LD = identity, V_SNP becomes diag of (SE * varSNP)^2 per SNP."""
    SE_SNP = np.array([[0.01, 0.02], [0.03, 0.04]])
    I_LD = np.eye(2)
    varSNP = np.array([0.4, 0.5])
    V = compute_v_snp_batch(SE_SNP, I_LD, varSNP)
    # Off-diagonals: I_LD[0,1] = 0 -> V[:,0,1] = 0
    np.testing.assert_allclose(V[:, 0, 1], 0.0)
    np.testing.assert_allclose(V[:, 1, 0], 0.0)
    # Diagonals: (SE * sqrt(I=1) * varSNP)^2
    expected_d0 = (SE_SNP[:, 0] * varSNP) ** 2
    expected_d1 = (SE_SNP[:, 1] * varSNP) ** 2
    np.testing.assert_allclose(V[:, 0, 0], expected_d0)
    np.testing.assert_allclose(V[:, 1, 1], expected_d1)


def test_v_snp_construction_with_intercept():
    """Verify off-diagonal formula: SE_x*SE_y*I_LD[x,y]*sqrt(I_xx)*sqrt(I_yy)*varSNP^2."""
    SE_SNP = np.array([[0.01, 0.02]])
    I_LD = np.array([[1.05, 0.20], [0.20, 1.10]])
    varSNP = np.array([0.4])
    V = compute_v_snp_batch(SE_SNP, I_LD, varSNP)
    expected_off = (
        SE_SNP[0, 0]
        * SE_SNP[0, 1]
        * I_LD[0, 1]
        * np.sqrt(I_LD[0, 0])
        * np.sqrt(I_LD[1, 1])
        * varSNP[0] ** 2
    )
    np.testing.assert_allclose(V[0, 0, 1], expected_off)
    np.testing.assert_allclose(V[0, 1, 0], expected_off)


def test_v_full_block_diagonal_structure():
    V_LD = np.diag([1e-4, 5e-5, 1.5e-4])
    V_SNP = np.array([[[1.0, 0.5], [0.5, 2.0]]])
    V_full = compute_v_full_batch(V_LD, varSNPSE2=2.5e-7, V_SNP_batch=V_SNP)
    assert V_full.shape == (1, 6, 6)
    # Top-left
    assert V_full[0, 0, 0] == pytest.approx(2.5e-7)
    # SNP block (rows/cols 1-2)
    np.testing.assert_array_equal(V_full[0, 1:3, 1:3], V_SNP[0])
    # LD block (rows/cols 3-5)
    np.testing.assert_array_equal(V_full[0, 3:, 3:], V_LD)
    # Off-diagonal blocks must be zero
    assert np.all(V_full[0, 0, 1:] == 0)
    assert np.all(V_full[0, 1:3, 3:] == 0)


def test_s_fullrun_layout():
    S_LD = np.array([[0.083, 0.034], [0.034, 0.070]])
    varSNP = np.array([0.4, 0.5])
    beta_SNP = np.array([[0.01, 0.005], [-0.02, 0.01]])
    S = compute_s_fullrun_batch(S_LD, varSNP, beta_SNP)
    assert S.shape == (2, 3, 3)
    np.testing.assert_allclose(S[:, 0, 0], varSNP)
    np.testing.assert_allclose(S[:, 0, 1], varSNP * beta_SNP[:, 0])
    np.testing.assert_allclose(S[:, 0, 2], varSNP * beta_SNP[:, 1])
    np.testing.assert_array_equal(S[0, 1:, 1:], S_LD)
    np.testing.assert_array_equal(S[1, 1:, 1:], S_LD)


def test_jacobian_via_finite_difference():
    """Numerical derivative of model_implied_vech matches model_jacobian."""
    rng = np.random.default_rng(42)
    N = 5
    theta = np.column_stack(
        [
            rng.uniform(0.5, 1.5, N),  # lambda
            rng.uniform(-0.05, 0.05, N),  # beta
            rng.uniform(0.01, 0.1, N),  # psi_1
            rng.uniform(0.01, 0.1, N),  # psi_2
            rng.uniform(0.05, 0.2, N),  # sigma_sq
        ]
    )
    varSNP = rng.uniform(0.1, 0.5, N)
    J = model_jacobian(theta, varSNP)
    h = 1e-6
    for j in range(5):
        bump = np.zeros_like(theta)
        bump[:, j] = h
        plus = model_implied_vech(theta + bump, varSNP)
        minus = model_implied_vech(theta - bump, varSNP)
        finite = (plus - minus) / (2 * h)
        np.testing.assert_allclose(J[:, :, j], finite, atol=1e-6, rtol=1e-4)


def test_fit_recovers_known_truth():
    """
    Build clean S_LD, S_Fullrun consistent with a known theta, fit, and
    verify recovery to high precision.
    """
    rng = np.random.default_rng(0)
    # True parameters
    lam_true = 0.85
    beta_true = 0.02
    psi_1_true = 0.04
    psi_2_true = 0.05
    sigma_sq_true = 0.06

    Lambda = np.array([1.0, lam_true])
    var_F_marginal = sigma_sq_true  # at beta=0 baseline
    # Trait-trait genetic covariance under truth (without SNP)
    S_LD = np.outer(Lambda, Lambda) * var_F_marginal + np.diag([psi_1_true, psi_2_true])

    N = 40
    varSNP = rng.uniform(0.1, 0.5, N)
    # beta_SNP for trait i should be lambda_i * beta_true at the truth
    beta_SNP = np.column_stack(
        [
            np.full(N, 1.0 * beta_true),
            np.full(N, lam_true * beta_true),
        ]
    )
    SE_SNP = np.full((N, 2), 0.01)
    I_LD = np.eye(2)

    # V_LD: small diagonal, well-conditioned
    V_LD = np.diag([1e-4, 5e-5, 1e-4])

    result = fit_common_factor_gwas(
        S_LD=S_LD,
        V_LD=V_LD,
        I_LD=I_LD,
        beta_SNP=beta_SNP,
        SE_SNP=SE_SNP,
        varSNP=varSNP,
        max_iter=50,
        tol=1e-12,
    )
    assert np.all(result.fail == 0), f"Failures: {result.fail}"
    np.testing.assert_allclose(result.est, beta_true, atol=1e-5)
    # SE should be finite & positive
    assert np.all(result.se_c > 0)
    assert np.all(np.isfinite(result.se_c))


# ---- R comparison test -----------------------------------------------------


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
def test_matches_lavaan_on_small_slice():
    """
    Run R's commonfactorGWAS and our kernel on the same well-conditioned
    synthetic data, and assert the per-SNP estimates and sandwich SEs match.
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    gsem = importr("GenomicSEM")

    rng = np.random.default_rng(42)
    N_SNPS = 30

    # Synthetic ground-truth common-factor model
    lam_true = 0.80
    psi_1_true = 0.04
    psi_2_true = 0.05
    sigma_sq_true = 0.06
    Lambda = np.array([1.0, lam_true])
    S_LD = np.outer(Lambda, Lambda) * sigma_sq_true + np.diag([psi_1_true, psi_2_true])

    # Reasonable LDSC sampling covariance: diagonal with small entries
    V_LD = np.diag([3e-5, 1e-5, 4e-5])
    I_LD = np.eye(2) * 1.05  # mild inflation

    # Per-SNP: random MAF, modest betas, realistic SE
    maf = rng.uniform(0.05, 0.5, N_SNPS)
    varSNP = 2.0 * maf * (1.0 - maf)
    # Generate per-SNP true betas around 0 with realistic spread
    beta_factor_true = rng.normal(0.0, 0.01, N_SNPS)  # SNP -> F1 effects
    # Implied per-trait OLS betas under the model: beta_trait_j = lambda_j * beta_F
    beta_SNP = np.column_stack([beta_factor_true, lam_true * beta_factor_true])
    # Add small observational noise to mimic LDSC noise on the SNP entries
    beta_SNP = beta_SNP + rng.normal(0.0, 0.001, beta_SNP.shape)
    SE_SNP = np.full((N_SNPS, 2), 0.01)

    # ---- Build R inputs ----
    base = importr("base")
    trait_names = ["t1", "t2"]
    r_S = base.matrix(ro.FloatVector(S_LD.flatten("F")), nrow=2)
    r_S.rownames = ro.StrVector(trait_names)
    r_S.colnames = ro.StrVector(trait_names)
    r_V = base.matrix(ro.FloatVector(V_LD.flatten("F")), nrow=3)
    r_I = base.matrix(ro.FloatVector(I_LD.flatten("F")), nrow=2)
    r_I.rownames = ro.StrVector(trait_names)
    r_I.colnames = ro.StrVector(trait_names)
    covstruc = ro.ListVector({"V": r_V, "S": r_S, "I": r_I})

    # SNPs dataframe in the format sumstats() emits
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

    r_result = gsem.commonfactorGWAS(
        covstruc=covstruc,
        SNPs=r_snps,
        parallel=False,
        toler=1e-50,  # relax singularity check; we want timing not failure
    )
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_result_df = ro.conversion.get_conversion().rpy2py(r_result)

    # ---- Run Python kernel ----
    py_result = fit_common_factor_gwas(
        S_LD=S_LD,
        V_LD=V_LD,
        I_LD=I_LD,
        beta_SNP=beta_SNP,
        SE_SNP=SE_SNP,
        varSNP=varSNP,
    )

    # ---- Compare ----
    # R returns one row per SNP in the same order as input; column 'est' is
    # the SNP -> F1 estimate, 'se_c' is the sandwich SE.
    r_est = np.asarray(r_result_df["est"], dtype=float)
    r_se_c = np.asarray(r_result_df["se_c"], dtype=float)

    # Restrict comparison to SNPs both sides flagged as ok
    ok = (py_result.fail == 0) & np.isfinite(r_est) & np.isfinite(r_se_c)
    assert ok.sum() >= 0.9 * N_SNPS, f"Too few comparable SNPs: {ok.sum()}/{N_SNPS}"

    # Empirically observed agreement on this fixture is much tighter than
    # lavaan's nominal optim.dx.tol=0.01 would suggest -- ~1e-9 absolute on
    # est and ~1e-7 on se_c -- so set tolerances accordingly. If lavaan ever
    # tightens its convergence criterion these can come down further.
    np.testing.assert_allclose(py_result.est[ok], r_est[ok], atol=1e-7, rtol=1e-5)
    np.testing.assert_allclose(py_result.se_c[ok], r_se_c[ok], atol=1e-5, rtol=1e-3)
