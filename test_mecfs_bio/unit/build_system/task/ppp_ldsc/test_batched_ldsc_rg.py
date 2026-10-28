"""Validate the batched cross-trait rg kernel against the repo's exact GenomicSEM port.

The batched kernel encodes per-protein filtering as zero weights over the shared SNP set, so its
gcov / h2_trait / h2_protein POINT estimates must be machine-identical to the exact per-pair
estimators (estimate_h2 + estimate_cov), and hence rg must match. The jackknife SE differs only
because block membership shifts by the few zero-weighted SNPs, so it is checked for finiteness /
plausibility rather than machine equality. Mirrors test_batched_ldsc_h2.
"""

import numpy as np
import pytest

from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_h2 import exact_h2_single
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_rg import (
    batched_rg,
    estimate_trait_context,
    exact_rg_single,
)

_N_BLOCKS = 20
_M = 100_000.0
_N = 10_000.0  # chi^2 filter threshold = max(0.001*N, 80) = 80


def _synthetic(
    seed: int, n_snps: int = 3000, k: int = 3
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    ld = rng.uniform(1.0, 50.0, size=n_snps)
    # Correlated trait/protein z: shared component + independent noise, scaled by sqrt(ld) so
    # chi^2 grows with LD as LDSC expects.
    shared = rng.standard_normal(n_snps) * np.sqrt(0.02 * ld)
    z_trait = (
        shared
        + rng.standard_normal(n_snps) * np.sqrt(0.02 * ld)
        + rng.standard_normal(n_snps)
    )
    z_protein = np.empty((n_snps, k))
    for j in range(k):
        z_protein[:, j] = (
            (0.5 + 0.2 * j) * shared
            + rng.standard_normal(n_snps) * np.sqrt(0.02 * ld)
            + rng.standard_normal(n_snps)
        )
    # A few missing variants per protein and in the trait.
    z_protein[rng.integers(0, n_snps, 8), 1] = np.nan
    z_trait[rng.integers(0, n_snps, 6)] = np.nan
    n_trait = np.full(n_snps, _N)
    return ld, z_trait, n_trait, z_protein


def test_batched_matches_exact_point_estimates() -> None:
    ld, z_trait, n_trait, z_protein = _synthetic(seed=0)
    k = z_protein.shape[1]
    n_protein = np.full(k, _N)
    trait_ctx = estimate_trait_context(
        z_trait=z_trait, n_trait=n_trait, ld=ld, m=_M, n_blocks=_N_BLOCKS
    )
    res = batched_rg(
        trait_ctx,
        z_protein=z_protein,
        ld=ld,
        n_protein=n_protein,
        m=_M,
        n_blocks=_N_BLOCKS,
    )

    for j in range(k):
        exact = exact_rg_single(
            z_trait, _N, z_protein[:, j], ld, _N, _M, n_blocks=_N_BLOCKS
        )
        assert np.isclose(res.h2_trait, exact.h2_trait, atol=1e-9)
        assert np.isclose(res.h2_protein[j], exact.h2_protein, atol=1e-9)
        assert np.isclose(res.gcov[j], exact.gcov, atol=1e-9)
        assert np.isclose(res.gcov_intercept[j], exact.gcov_intercept, atol=1e-9)
        assert np.isclose(res.rg[j], exact.rg, atol=1e-9)
        # SE: finite, positive, plausible magnitude.
        assert np.isfinite(res.rg_se[j]) and res.rg_se[j] > 0


def test_batched_exclude_mask_matches_exact_with_exclusion() -> None:
    ld, z_trait, n_trait, z_protein = _synthetic(seed=1)
    k = z_protein.shape[1]
    n_protein = np.full(k, _N)
    exclude = np.zeros(z_protein.shape, dtype=bool)
    exclude[1000:1200, 0] = True  # a contiguous "cis" block for protein 0 only

    trait_ctx = estimate_trait_context(
        z_trait=z_trait, n_trait=n_trait, ld=ld, m=_M, n_blocks=_N_BLOCKS
    )
    res = batched_rg(
        trait_ctx,
        z_protein=z_protein,
        ld=ld,
        n_protein=n_protein,
        m=_M,
        n_blocks=_N_BLOCKS,
        exclude=exclude,
    )
    exact0 = exact_rg_single(
        z_trait,
        _N,
        z_protein[:, 0],
        ld,
        _N,
        _M,
        n_blocks=_N_BLOCKS,
        exclude=exclude[:, 0],
    )
    assert np.isclose(res.h2_protein[0], exact0.h2_protein, atol=1e-9)
    assert np.isclose(res.gcov[0], exact0.gcov, atol=1e-9)
    assert np.isclose(res.rg[0], exact0.rg, atol=1e-9)
    # Excluded variants reduce the covariance SNP count for protein 0 only.
    assert res.n_snps[0] < res.n_snps[2]


def test_rg_reports_reasonable_diagnostics() -> None:
    ld, z_trait, n_trait, z_protein = _synthetic(seed=2)
    k = z_protein.shape[1]
    n_protein = np.full(k, _N)
    trait_ctx = estimate_trait_context(
        z_trait=z_trait, n_trait=n_trait, ld=ld, m=_M, n_blocks=_N_BLOCKS
    )
    res = batched_rg(
        trait_ctx,
        z_protein=z_protein,
        ld=ld,
        n_protein=n_protein,
        m=_M,
        n_blocks=_N_BLOCKS,
    )
    assert np.all(res.n_snps > 2900)
    assert np.all(np.isfinite(res.rg))
    assert np.all(np.isfinite(res.rg_p))
    # The proteins share a positive component with the trait, so rg is positive and bounded.
    assert np.all(res.rg > 0)
    assert np.all(np.abs(res.rg) < 1.5)


def test_trait_context_reports_h2_se_and_intercept_matching_exact() -> None:
    rng = np.random.default_rng(3)
    s = 3000
    ld = rng.uniform(1.0, 50.0, size=s)
    z = rng.standard_normal(s) * np.sqrt(0.02 * ld) + rng.standard_normal(s)
    ctx = estimate_trait_context(
        z_trait=z, n_trait=np.full(s, _N), ld=ld, m=_M, n_blocks=_N_BLOCKS
    )
    exact = exact_h2_single(z**2, ld, _N, _M, n_blocks=_N_BLOCKS)
    assert ctx.h2 == pytest.approx(exact.h2, rel=1e-8)
    assert ctx.intercept == pytest.approx(exact.intercept, rel=1e-8)
    # Block membership differs slightly between the shared-block and exact jackknives.
    assert ctx.h2_se == pytest.approx(exact.h2_se, rel=0.05)
    assert ctx.intercept_se > 0
    assert ctx.n_snps == int(ctx.keep.sum())
    assert ctx.mean_chi2 > 1.0
    assert ctx.lambda_gc > 1.0


def test_rg_reports_finite_gcov_and_intercept_ses() -> None:
    ld, z_trait, n_trait, z_protein = _synthetic(seed=4)
    n_protein = np.full(z_protein.shape[1], _N)
    ctx = estimate_trait_context(
        z_trait=z_trait, n_trait=n_trait, ld=ld, m=_M, n_blocks=_N_BLOCKS
    )
    res = batched_rg(
        ctx, z_protein=z_protein, ld=ld, n_protein=n_protein, m=_M, n_blocks=_N_BLOCKS
    )
    assert np.all(np.isfinite(res.gcov_se)) and np.all(res.gcov_se > 0)
    assert np.all(np.isfinite(res.gcov_intercept_se))
    assert np.all(res.gcov_intercept_se > 0)


def _per_snp_n_traits(
    seed: int, n_snps: int = 3000
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Three correlated traits with per-SNP N and about 5% missing variants per column."""
    rng = np.random.default_rng(seed)
    ld = rng.uniform(1.0, 50.0, size=n_snps)
    shared = rng.standard_normal(n_snps) * np.sqrt(0.02 * ld)
    z = np.empty((n_snps, 3))
    for j in range(3):
        z[:, j] = (
            shared
            + rng.standard_normal(n_snps) * np.sqrt(0.02 * ld)
            + rng.standard_normal(n_snps)
        )
    n = rng.uniform(8_000.0, 12_000.0, size=(n_snps, 3))
    missing = rng.random((n_snps, 3)) < 0.05
    z[missing] = np.nan
    n[missing] = np.nan
    return ld, z, n


def test_per_snp_n_second_side_matches_its_own_trait_context() -> None:
    ld, z, n = _per_snp_n_traits(seed=5)
    ctx0 = estimate_trait_context(
        z_trait=z[:, 0], n_trait=n[:, 0], ld=ld, m=_M, n_blocks=_N_BLOCKS
    )
    res = batched_rg(
        ctx0, z_protein=z[:, 1:], ld=ld, n_protein=n[:, 1:], m=_M, n_blocks=_N_BLOCKS
    )
    for j in (1, 2):
        own = estimate_trait_context(
            z_trait=z[:, j], n_trait=n[:, j], ld=ld, m=_M, n_blocks=_N_BLOCKS
        )
        assert res.h2_protein[j - 1] == pytest.approx(own.h2, rel=1e-10)
        assert res.n_bar_protein[j - 1] == pytest.approx(own.n_bar, rel=1e-12)


def test_constant_matrix_n_equals_vector_n() -> None:
    ld, z_trait, n_trait, z_protein = _synthetic(seed=6)
    ctx = estimate_trait_context(
        z_trait=z_trait, n_trait=n_trait, ld=ld, m=_M, n_blocks=_N_BLOCKS
    )
    k = z_protein.shape[1]
    vec = batched_rg(
        ctx,
        z_protein=z_protein,
        ld=ld,
        n_protein=np.full(k, _N),
        m=_M,
        n_blocks=_N_BLOCKS,
    )
    mat = batched_rg(
        ctx,
        z_protein=z_protein,
        ld=ld,
        n_protein=np.where(np.isfinite(z_protein), _N, np.nan),
        m=_M,
        n_blocks=_N_BLOCKS,
    )
    np.testing.assert_allclose(mat.rg, vec.rg, rtol=1e-12)
    np.testing.assert_allclose(mat.rg_se, vec.rg_se, rtol=1e-12)
