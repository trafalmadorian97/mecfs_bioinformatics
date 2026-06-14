"""
Pure-Python tests for the polars LD-score munge port (`_genomic_sem_munge`):
the odds-ratio guard and the N / NaN-handling regression tests. The numerical
comparison against GenomicSEM::munge lives in
test_genomic_sem_munge_r_comparison.py.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_munge import (
    munge_sumstats,
)


def test_odds_ratio_effect_raises():
    """
    An effect column that looks like an odds ratio (median rounds to 1) is a
    mis-specified input here: this pipeline always feeds a log-scale beta
    (gwaslab keeps ORs in a separate column). Unlike GenomicSEM, which silently
    log-transforms, munge_sumstats raises rather than mutate the data.
    """
    rng = np.random.default_rng(11)
    n_snps = 300
    snps = [f"rs{i}" for i in range(n_snps)]
    bases = np.array(["A", "C", "G", "T"])
    a1 = rng.choice(bases, n_snps)
    a2 = np.array([rng.choice([b for b in bases if b != a]) for a in a1])
    odds_ratio = np.exp(rng.normal(0, 0.1, n_snps))  # median ~ 1

    sumstats = pl.DataFrame(
        {
            "SNP": snps,
            "A1": a1,
            "A2": a2,
            "effect": odds_ratio,
            "SE": np.abs(rng.normal(0.05, 0.01, n_snps)),
            "P": rng.uniform(1e-6, 1.0, n_snps),
            "N": 9999.0,
        }
    )
    ref = pl.DataFrame({"SNP": snps, "A1": a1, "A2": a2})

    with pytest.raises(ValueError, match="looks like an odds ratio"):
        munge_sumstats(sumstats, ref, n=9999.0)


def _simple_matching_inputs() -> tuple[pl.DataFrame, pl.DataFrame]:
    """Three SNPs whose file alleles match the reference, with a per-SNP N."""
    snps = ["rs0", "rs1", "rs2"]
    sumstats = pl.DataFrame(
        {
            "SNP": snps,
            "A1": ["A", "C", "G"],
            "A2": ["G", "T", "A"],
            "effect": [0.1, -0.2, 0.05],
            "SE": [0.05, 0.05, 0.05],
            "P": [0.01, 0.2, 0.5],
            "N": [11000.0, 12000.0, 13000.0],
        }
    )
    ref = pl.DataFrame({"SNP": snps, "A1": ["A", "C", "G"], "A2": ["G", "T", "A"]})
    return sumstats, ref


def test_nan_n_keeps_file_n_column():
    """
    A NaN sample size means "not provided" (mirroring R's `!is.na(N)`), so the
    file's per-SNP N column must survive rather than be clobbered with NaN.
    Regression test for the bug where `if n is not None` overrode N with NaN.
    """
    sumstats, ref = _simple_matching_inputs()
    out = munge_sumstats(sumstats, ref, n=float("nan")).sort("SNP")
    np.testing.assert_array_equal(
        out["N"].to_numpy(), np.array([11000.0, 12000.0, 13000.0])
    )


def test_scalar_n_overrides_file_n_column():
    """A real scalar N overrides the file's N column, as GenomicSEM does."""
    sumstats, ref = _simple_matching_inputs()
    out = munge_sumstats(sumstats, ref, n=50000.0).sort("SNP")
    np.testing.assert_array_equal(out["N"].to_numpy(), np.full(3, 50000.0))


def test_nan_p_or_effect_rows_are_dropped():
    """
    R's is.na() drops both NA and NaN, so a NaN P or NaN effect row must be
    dropped (polars is_not_null alone would keep NaN). Regression test.
    """
    sumstats, ref = _simple_matching_inputs()
    sumstats = sumstats.with_columns(
        pl.Series("P", [float("nan"), 0.2, 0.5]),
        pl.Series("effect", [0.1, float("nan"), 0.05]),
    )
    out = munge_sumstats(sumstats, ref, n=50000.0).sort("SNP")
    # rs0 (NaN P) and rs1 (NaN effect) drop; only rs2 survives.
    assert out["SNP"].to_list() == ["rs2"]
