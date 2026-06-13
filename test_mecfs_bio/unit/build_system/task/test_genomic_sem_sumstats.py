"""
Tests for the polars sumstats port (`_genomic_sem_sumstats`).

Headline comparison against GenomicSEM::sumstats on synthetic data with an
se.logit trait (carrying its own MAF) and an OLS trait (using the reference
MAF), covering allele flips/mismatches, MAF folding, and listwise merge.
"""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_sumstats import (
    SumstatsTrait,
    run_sumstats,
)


def _have_rpy2_genomic_sem() -> bool:
    try:
        import rpy2.robjects as ro  # noqa: F401
        from rpy2.robjects.packages import importr

        importr("GenomicSEM")
        return True
    except Exception:
        return False


def _make_reference(rng, n_snps: int) -> tuple[pl.DataFrame, np.ndarray, np.ndarray]:
    bases = np.array(["A", "C", "G", "T"])
    a1 = rng.choice(bases, n_snps)
    a2 = np.array([rng.choice([b for b in bases if b != a]) for a in a1])
    ref = pl.DataFrame(
        {
            "SNP": [f"rs{i}" for i in range(n_snps)],
            "CHR": rng.integers(1, 23, n_snps),
            "BP": rng.integers(1, 100_000_000, n_snps),
            "MAF": rng.uniform(0.005, 0.5, n_snps),
            "A1": a1,
            "A2": a2,
        }
    )
    return ref, a1, a2


def _trait_alleles(rng, a1_ref, a2_ref):
    """Mostly-matching file alleles with some swaps and mismatches."""
    bases = np.array(["A", "C", "G", "T"])
    n = len(a1_ref)
    kind = rng.choice(["match", "swap", "mismatch"], n, p=[0.7, 0.2, 0.1])
    a1 = a1_ref.copy()
    a2 = a2_ref.copy()
    for i, k in enumerate(kind):
        if k == "swap":
            a1[i], a2[i] = a2_ref[i], a1_ref[i]
        elif k == "mismatch":
            others = [b for b in bases if b not in (a1_ref[i], a2_ref[i])]
            a1[i] = others[0]
            a2[i] = others[1]
    return a1, a2


def _run_r_sumstats(tmp_path, ref, trait_files, trait_names, ols, se_logit, ns):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    gsem = importr("GenomicSEM")
    base = importr("base")

    ref_path = tmp_path / "ref.txt"
    ref.write_csv(ref_path, separator="\t")
    paths = []
    for name, df in zip(trait_names, trait_files):
        p = tmp_path / f"{name}.txt"
        df.write_csv(p, separator="\t")
        paths.append(str(p))

    cwd = str(base.getwd()[0])
    try:
        base.setwd(str(tmp_path))
        r_df = gsem.sumstats(
            files=ro.StrVector(paths),
            ref=str(ref_path),
            trait_names=ro.StrVector(trait_names),
            se_logit=ro.BoolVector(se_logit),
            OLS=ro.BoolVector(ols),
            linprob=ro.BoolVector([False] * len(trait_names)),
            N=ro.FloatVector(ns),
            info_filter=0.6,
            maf_filter=0.01,
            ambig=False,
        )
    finally:
        base.setwd(cwd)
    with localconverter(ro.default_converter + pandas2ri.converter):
        return pl.from_pandas(ro.conversion.get_conversion().rpy2py(r_df))


@pytest.mark.skipif(
    not _have_rpy2_genomic_sem(), reason="rpy2/GenomicSEM not available"
)
def test_matches_genomic_sem_sumstats(tmp_path):
    rng = np.random.default_rng(5)
    n_snps = 500
    ref, a1_ref, a2_ref = _make_reference(rng, n_snps)

    # Trait 0: se.logit (logistic betas + logistic SE), carries its own MAF.
    a1_0, a2_0 = _trait_alleles(rng, a1_ref, a2_ref)
    t0 = pl.DataFrame(
        {
            "SNP": ref["SNP"],
            "A1": a1_0,
            "A2": a2_0,
            "effect": rng.normal(0, 0.1, n_snps),  # log OR
            "SE": np.abs(rng.normal(0.05, 0.01, n_snps)),
            "P": rng.uniform(1e-6, 1.0, n_snps),
            "N": 80000.0,
            "MAF": rng.uniform(0.02, 0.98, n_snps),  # own MAF, some > 0.5
        }
    )
    # Trait 1: OLS (continuous), no MAF column -> uses reference MAF.
    a1_1, a2_1 = _trait_alleles(rng, a1_ref, a2_ref)
    p1 = rng.uniform(1e-6, 1.0, n_snps)
    # P == 1 -> Z == 0 -> standardized effect 0 -> se NaN; both must drop these.
    p1[:10] = 1.0
    t1 = pl.DataFrame(
        {
            "SNP": ref["SNP"],
            "A1": a1_1,
            "A2": a2_1,
            "effect": rng.normal(0, 0.05, n_snps),
            "SE": np.abs(rng.normal(0.02, 0.005, n_snps)),
            "P": p1,
            "N": 120000.0,
        }
    )

    trait_names = ["T0", "T1"]
    r_out = _run_r_sumstats(
        tmp_path,
        ref,
        [t0, t1],
        trait_names,
        ols=[False, True],
        se_logit=[True, False],
        ns=[80000.0, 120000.0],
    ).sort("SNP")

    py_out = run_sumstats(
        [
            SumstatsTrait(df=t0, name="T0", n=80000.0, se_logit=True),
            SumstatsTrait(df=t1, name="T1", n=120000.0, ols=True),
        ],
        ref,
        maf_filter=0.01,
        info_filter=0.6,
    ).sort("SNP")

    assert r_out["SNP"].to_list() == py_out["SNP"].to_list()
    for col in ["beta.T0", "se.T0", "beta.T1", "se.T1", "MAF"]:
        np.testing.assert_allclose(
            py_out[col].to_numpy(), r_out[col].to_numpy(), rtol=1e-6, atol=1e-8
        )
    assert py_out["A1"].to_list() == r_out["A1"].to_list()
    assert py_out["A2"].to_list() == r_out["A2"].to_list()
