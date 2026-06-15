"""
Unit tests for the polars sumstats standardisation (`run_sumstats`). Every trait
is standardised on one linear scale (`beta = Z/sqrt(N*varSNP)`, `se = |beta/Z|`);
there is no binary/logistic special-casing. The odds-ratio guard is retained
because the linear transform relies on sign(effect).
"""

from __future__ import annotations

import polars as pl
import pytest
from scipy.stats import norm

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_sumstats import (
    SumstatsTrait,
    run_sumstats,
)


def _ref(snps, a1, a2, maf) -> pl.DataFrame:
    return pl.DataFrame({"SNP": snps, "A1": a1, "A2": a2, "MAF": maf})


def test_linear_standardization():
    """beta = Z/sqrt(N*varSNP) and se = 1/sqrt(N*varSNP), for any trait."""
    ref = _ref(["rs1"], ["A"], ["G"], [0.3])
    df = pl.DataFrame(
        {
            "SNP": ["rs1"],
            "A1": ["A"],
            "A2": ["G"],
            "effect": [0.1],  # log-scale beta; only its sign is used
            "SE": [0.04],
            "P": [0.05],
            "N": [10000.0],
        }
    )
    out = run_sumstats([SumstatsTrait(df=df, name="t", n=None)], ref)

    z = float(norm.isf(0.05 / 2.0))
    varsnp = 2 * 0.3 * (1 - 0.3)
    expected_beta = z / (10000.0 * varsnp) ** 0.5
    expected_se = 1.0 / (10000.0 * varsnp) ** 0.5

    assert out["beta.t"][0] == pytest.approx(expected_beta, rel=1e-9)
    assert out["se.t"][0] == pytest.approx(expected_se, rel=1e-9)


def test_negative_effect_flips_z_sign():
    """A negative effect gives a negative standardised beta (sign carried through)."""
    ref = _ref(["rs1"], ["A"], ["G"], [0.3])
    df = pl.DataFrame(
        {
            "SNP": ["rs1"],
            "A1": ["A"],
            "A2": ["G"],
            "effect": [-0.2],
            "SE": [0.04],
            "P": [0.05],
            "N": [10000.0],
        }
    )
    out = run_sumstats([SumstatsTrait(df=df, name="t", n=None)], ref)
    assert out["beta.t"][0] < 0


def test_odds_ratio_effect_raises():
    """An effect column whose median rounds to 1 looks like an OR -> raise."""
    ref = _ref(["rs1", "rs2", "rs3"], ["A", "C", "G"], ["G", "T", "A"], [0.3, 0.3, 0.3])
    df = pl.DataFrame(
        {
            "SNP": ["rs1", "rs2", "rs3"],
            "A1": ["A", "C", "G"],
            "A2": ["G", "T", "A"],
            "effect": [0.95, 1.0, 1.05],  # median rounds to 1
            "SE": [0.04, 0.04, 0.04],
            "P": [0.05, 0.2, 0.5],
            "N": [10000.0, 10000.0, 10000.0],
        }
    )
    with pytest.raises(ValueError, match="odds ratio"):
        run_sumstats([SumstatsTrait(df=df, name="t", n=None)], ref)
