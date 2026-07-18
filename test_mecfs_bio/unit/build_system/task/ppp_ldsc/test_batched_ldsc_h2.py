"""Validate the batched heritability kernel against the repo's exact GenomicSEM port.

The batched kernel encodes per-protein filtering as zero weights over the shared SNP set,
so its h2 POINT estimate must be machine-identical to the exact per-protein estimator, and
its jackknife SE must agree closely (differing only because block membership shifts by the
few zero-weighted SNPs). This mirrors experiments/claude/ppp_ldsc/batched_vs_exact_h2_probe.
"""

import numpy as np

from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_h2 import (
    batched_h2,
    exact_h2_single,
)

_N_BLOCKS = 20
_M = 100_000.0
_N = 10_000.0  # chi^2 filter threshold = max(0.001*N, 80) = 80


def _synthetic(seed: int, n_snps: int = 3000, k: int = 3):
    rng = np.random.default_rng(seed)
    ld = rng.uniform(1.0, 50.0, size=n_snps)
    chi2 = 1.0 + 0.02 * ld[:, None] + rng.exponential(1.0, size=(n_snps, k))
    # A few variants above the chi^2 filter, and a few missing, per protein.
    chi2[rng.integers(0, n_snps, 5), 0] = 200.0
    chi2[rng.integers(0, n_snps, 8), 1] = np.nan
    return ld, chi2


def test_batched_matches_exact_point_estimate_and_close_se():
    ld, chi2 = _synthetic(seed=0)
    n = np.full(chi2.shape[1], _N)
    batched = batched_h2(chi2, ld, n, _M, n_blocks=_N_BLOCKS)

    for j in range(chi2.shape[1]):
        exact = exact_h2_single(chi2[:, j], ld, _N, _M, n_blocks=_N_BLOCKS)
        # Point estimate + intercept: machine-identical (zero-weight == dropped).
        assert np.isclose(batched.h2[j], exact.h2, atol=1e-9)
        assert np.isclose(batched.intercept[j], exact.intercept, atol=1e-9)
        # Jackknife SE: close (shared-block vs kept-block membership).
        assert np.isclose(batched.h2_se[j], exact.h2_se, rtol=0.1)


def test_batched_exclude_mask_matches_exact_with_exclusion():
    ld, chi2 = _synthetic(seed=1)
    n = np.full(chi2.shape[1], _N)
    # Exclude a contiguous "cis" block for protein 0 only.
    exclude = np.zeros(chi2.shape, dtype=bool)
    exclude[1000:1200, 0] = True

    batched = batched_h2(chi2, ld, n, _M, n_blocks=_N_BLOCKS, exclude=exclude)
    exact0 = exact_h2_single(
        chi2[:, 0], ld, _N, _M, n_blocks=_N_BLOCKS, exclude=exclude[:, 0]
    )
    assert np.isclose(batched.h2[0], exact0.h2, atol=1e-9)
    # The excluded variants are reflected in the kept-SNP count.
    assert batched.n_snps[0] < batched.n_snps[1]


def test_batched_reports_reasonable_diagnostics():
    ld, chi2 = _synthetic(seed=2)
    n = np.full(chi2.shape[1], _N)
    res = batched_h2(chi2, ld, n, _M, n_blocks=_N_BLOCKS)
    assert np.all(res.n_snps > 2900)  # only a handful filtered
    assert np.all(res.mean_chi2 > 1.0)
    assert np.all(np.isfinite(res.lambda_gc))
