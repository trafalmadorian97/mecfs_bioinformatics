"""
Tests for the polars LD-score munge port (`_genomic_sem_munge`).

Headline comparison against GenomicSEM::munge on synthetic summary statistics
that exercise allele flips, allele mismatches, MAF folding, and INFO/MAF
filtering, plus a focused odds-ratio-detection test.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_munge import (
    munge_sumstats,
)


def _have_rpy2_genomic_sem() -> bool:
    try:
        import rpy2.robjects as ro  # noqa: F401
        from rpy2.robjects.packages import importr

        importr("GenomicSEM")
        return True
    except Exception:
        return False


def _make_inputs(seed: int) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Build a synthetic sumstats frame and HapMap3-style reference with a mix of
    matching, allele-swapped, mismatched, non-ACGT, and low-MAF/low-INFO SNPs.
    """
    rng = np.random.default_rng(seed)
    n_snps = 400
    snps = [f"rs{i}" for i in range(n_snps)]
    # Reference alleles.
    bases = np.array(["A", "C", "G", "T"])
    a1_ref = rng.choice(bases, n_snps)
    a2_ref = np.array([rng.choice([b for b in bases if b != a]) for a in a1_ref])

    # File alleles: most match, some swapped, some mismatched, some non-ACGT.
    kind = rng.choice(
        ["match", "swap", "mismatch", "bad"], n_snps, p=[0.6, 0.25, 0.1, 0.05]
    )
    a1_file = a1_ref.copy()
    a2_file = a2_ref.copy()
    for i, k in enumerate(kind):
        if k == "swap":
            a1_file[i], a2_file[i] = a2_ref[i], a1_ref[i]
        elif k == "mismatch":
            others = [b for b in bases if b not in (a1_ref[i], a2_ref[i])]
            a1_file[i] = others[0]
            a2_file[i] = others[1] if len(others) > 1 else others[0]
        elif k == "bad":
            a1_file[i] = "N"

    effect = rng.normal(0, 0.1, n_snps)
    p = rng.uniform(1e-6, 1.0, n_snps)
    maf = rng.uniform(0.001, 0.99, n_snps)  # some below 0.01, some above 0.5
    info = rng.uniform(0.5, 1.0, n_snps)  # some below 0.9

    sumstats = pl.DataFrame(
        {
            "SNP": snps,
            "A1": a1_file,
            "A2": a2_file,
            "effect": effect,
            "SE": np.abs(rng.normal(0.05, 0.01, n_snps)),
            "P": p,
            "N": 12345.0,
            "MAF": maf,
            "INFO": info,
        }
    )
    ref = pl.DataFrame({"SNP": snps, "A1": a1_ref, "A2": a2_ref})
    return sumstats, ref


def _run_r_munge(tmp_path: Path, sumstats: pl.DataFrame, ref: pl.DataFrame, n: float):
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    gsem = importr("GenomicSEM")
    base = importr("base")

    in_path = tmp_path / "trait.sumstats.txt"
    ref_path = tmp_path / "ref.snplist"
    sumstats.write_csv(in_path, separator="\t")
    ref.write_csv(ref_path, separator="\t")

    cwd = str(base.getwd()[0])
    try:
        base.setwd(str(tmp_path))
        gsem.munge(
            files=ro.StrVector([str(in_path)]),
            hm3=str(ref_path),
            trait_names=ro.StrVector(["trait"]),
            N=ro.FloatVector([n]),
            info_filter=0.9,
            maf_filter=0.01,
        )
    finally:
        base.setwd(cwd)
    return pl.read_csv(tmp_path / "trait.sumstats.gz", separator="\t")


@pytest.mark.skipif(
    not _have_rpy2_genomic_sem(), reason="rpy2/GenomicSEM not available"
)
def test_matches_genomic_sem_munge(tmp_path):
    sumstats, ref = _make_inputs(seed=3)
    n = 12345.0

    r_out = _run_r_munge(tmp_path, sumstats, ref, n)
    py_out = munge_sumstats(sumstats, ref, n=n, info_filter=0.9, maf_filter=0.01)

    r_sorted = r_out.sort("SNP")
    py_sorted = py_out.sort("SNP")

    assert r_sorted["SNP"].to_list() == py_sorted["SNP"].to_list()
    assert r_sorted["A1"].to_list() == py_sorted["A1"].to_list()
    assert r_sorted["A2"].to_list() == py_sorted["A2"].to_list()
    np.testing.assert_allclose(
        py_sorted["N"].to_numpy(), r_sorted["N"].to_numpy(), rtol=0, atol=0
    )
    np.testing.assert_allclose(
        py_sorted["Z"].to_numpy(), r_sorted["Z"].to_numpy(), rtol=1e-6, atol=1e-8
    )


@pytest.mark.skipif(
    not _have_rpy2_genomic_sem(), reason="rpy2/GenomicSEM not available"
)
def test_matches_genomic_sem_munge_odds_ratio(tmp_path):
    """When the effect column is an odds ratio (median ~ 1), both take log."""
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

    r_out = _run_r_munge(tmp_path, sumstats, ref, 9999.0).sort("SNP")
    py_out = munge_sumstats(sumstats, ref, n=9999.0).sort("SNP")

    assert r_out["SNP"].to_list() == py_out["SNP"].to_list()
    np.testing.assert_allclose(
        py_out["Z"].to_numpy(), r_out["Z"].to_numpy(), rtol=1e-6, atol=1e-8
    )


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
