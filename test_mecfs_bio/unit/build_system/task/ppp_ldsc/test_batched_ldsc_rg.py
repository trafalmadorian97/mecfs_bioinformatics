"""Validate the batched cross-trait rg kernel against the repo's exact GenomicSEM port.

The batched kernel encodes per-protein filtering as zero weights over the shared SNP set, so its
gcov / h2_trait / h2_protein POINT estimates must be machine-identical to the exact per-pair
estimators (estimate_h2 + estimate_cov), and hence rg must match. The jackknife SE differs only
because block membership shifts by the few zero-weighted SNPs, so it is checked for finiteness /
plausibility rather than machine equality. Mirrors test_batched_ldsc_h2.
"""

import numpy as np

from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_rg import (
    batched_rg,
    estimate_trait_context,
    exact_rg_single,
)

_N_BLOCKS = 20
_M = 100_000.0
_N = 10_000.0  # chi^2 filter threshold = max(0.001*N, 80) = 80


def _synthetic(seed: int, n_snps: int = 3000, k: int = 3):
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


def test_batched_matches_exact_point_estimates():
    ld, z_trait, n_trait, z_protein = _synthetic(seed=0)
    k = z_protein.shape[1]
    n_protein = np.full(k, _N)
    trait_ctx = estimate_trait_context(z_trait, n_trait, ld, _M, n_blocks=_N_BLOCKS)
    res = batched_rg(trait_ctx, z_protein, ld, n_protein, _M, n_blocks=_N_BLOCKS)

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


def test_batched_exclude_mask_matches_exact_with_exclusion():
    ld, z_trait, n_trait, z_protein = _synthetic(seed=1)
    k = z_protein.shape[1]
    n_protein = np.full(k, _N)
    exclude = np.zeros(z_protein.shape, dtype=bool)
    exclude[1000:1200, 0] = True  # a contiguous "cis" block for protein 0 only

    trait_ctx = estimate_trait_context(z_trait, n_trait, ld, _M, n_blocks=_N_BLOCKS)
    res = batched_rg(
        trait_ctx, z_protein, ld, n_protein, _M, n_blocks=_N_BLOCKS, exclude=exclude
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


def test_rg_reports_reasonable_diagnostics():
    ld, z_trait, n_trait, z_protein = _synthetic(seed=2)
    k = z_protein.shape[1]
    n_protein = np.full(k, _N)
    trait_ctx = estimate_trait_context(z_trait, n_trait, ld, _M, n_blocks=_N_BLOCKS)
    res = batched_rg(trait_ctx, z_protein, ld, n_protein, _M, n_blocks=_N_BLOCKS)
    assert np.all(res.n_snps > 2900)
    assert np.all(np.isfinite(res.rg))
    assert np.all(np.isfinite(res.rg_p))
    # The proteins share a positive component with the trait, so rg is positive and bounded.
    assert np.all(res.rg > 0)
    assert np.all(np.abs(res.rg) < 1.5)
