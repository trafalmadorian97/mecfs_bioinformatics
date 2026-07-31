"""
Pure-numpy unit tests for the LD-score regression building blocks
(`genomic_sem_ldsc`). The headline comparison against GenomicSEM::ldsc lives in
test_genomic_sem_ldsc_r_comparison.py.
"""

from __future__ import annotations

import numpy as np
import pytest

from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_ldsc import (
    _liability_conversion_factor,
    _regress_jackknife,
    block_bounds,
)


def test_block_bounds_partition_is_contiguous_and_complete():
    n, n_blocks = 1000, 50
    bounds = block_bounds(n, n_blocks)
    assert len(bounds) == n_blocks
    assert bounds[0][0] == 0
    assert bounds[-1][1] == n
    # contiguous, non-overlapping cover of [0, n)
    for (a, b), (c, _d) in zip(bounds, bounds[1:]):
        assert a < b == c
    covered = sum(b - a for a, b in bounds)
    assert covered == n


def test_regress_jackknife_recovers_known_line():
    rng = np.random.default_rng(0)
    n = 2000
    ld = rng.uniform(1, 100, n)
    true_slope, true_intercept = 3e-4, 1.02
    chi = true_intercept + true_slope * ld + rng.normal(0, 0.05, n)
    design = np.column_stack([ld, np.ones(n)])
    reg, pseudo = _regress_jackknife(design, chi, n_blocks=100)
    assert reg[0] == pytest.approx(true_slope, abs=5e-5)
    assert reg[1] == pytest.approx(true_intercept, abs=5e-2)
    assert pseudo.shape == (100, 2)
    # Jackknife mean of pseudo-values approximates the full estimate.
    np.testing.assert_allclose(pseudo.mean(axis=0), reg, rtol=0.2)


def test_liability_conversion_factor():
    assert _liability_conversion_factor(None, None) == 1.0
    assert _liability_conversion_factor(0.5, float("nan")) == 1.0
    # Known value for a binary trait.
    from scipy.stats import norm

    samp, pop = 0.3, 0.01
    z = norm.pdf(norm.ppf(1 - pop))
    expected = (pop**2 * (1 - pop) ** 2) / (samp * (1 - samp) * z**2)
    assert _liability_conversion_factor(samp, pop) == pytest.approx(expected)
