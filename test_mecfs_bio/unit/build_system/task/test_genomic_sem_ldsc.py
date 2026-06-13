"""
Tests for the pure-Python LD-score regression (`_genomic_sem_ldsc`).

Pure-numpy unit tests for the building blocks, plus a headline comparison
against GenomicSEM::ldsc on synthetic LD scores + munged sumstats.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_ldsc import (
    _block_bounds,
    _liability_conversion_factor,
    _regress_jackknife,
    run_ldsc,
)

# ---- pure-numpy unit tests --------------------------------------------------


def test_block_bounds_partition_is_contiguous_and_complete():
    n, n_blocks = 1000, 50
    bounds = _block_bounds(n, n_blocks)
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


# ---- R comparison test ------------------------------------------------------


def _have_rpy2_genomic_sem() -> bool:
    try:
        import rpy2.robjects as ro  # noqa: F401
        from rpy2.robjects.packages import importr

        importr("GenomicSEM")
        return True
    except Exception:
        return False


def _write_synthetic_inputs(
    tmp_dir: Path, *, n_chrom: int, snps_per_chrom: int, seed: int
) -> tuple[list[Path], Path, tuple[float, float], tuple[float, float]]:
    """
    Generate synthetic LD scores (<chr>.l2.ldscore.gz, <chr>.l2.M_5_50) and two
    munged sumstats files with a genetic signal and a positive cross-trait
    correlation. Returns (munged_paths, ld_dir, sample_prev, population_prev).
    """
    rng = np.random.default_rng(seed)
    ld_dir = tmp_dir / "ld"
    ld_dir.mkdir(parents=True, exist_ok=True)

    chrom_arr, snp_arr, bp_arr, l2_arr = [], [], [], []
    snp_counter = 0
    # Realistic reference SNP count so N*h2*L2/M stays O(1) (mild inflation);
    # otherwise the chi^2 > 80 filter would gut the SNP set.
    m_per_chrom = 500_000
    for chrom in range(1, n_chrom + 1):
        bp = np.sort(
            rng.choice(np.arange(1, 50_000_000), snps_per_chrom, replace=False)
        )
        l2 = rng.uniform(1.0, 80.0, snps_per_chrom)
        snps = [f"rs{snp_counter + i}" for i in range(snps_per_chrom)]
        snp_counter += snps_per_chrom
        ld_df = pd.DataFrame(
            {
                "CHR": chrom,
                "SNP": snps,
                "BP": bp,
                "CM": 0.0,
                "MAF": rng.uniform(0.05, 0.5, snps_per_chrom),
                "L2": l2,
            }
        )
        ld_df.to_csv(
            ld_dir / f"{chrom}.l2.ldscore.gz", sep="\t", index=False, compression="gzip"
        )
        (ld_dir / f"{chrom}.l2.M_5_50").write_text(f"{m_per_chrom}\n")
        chrom_arr.append(np.full(snps_per_chrom, chrom))
        snp_arr.extend(snps)
        bp_arr.append(bp)
        l2_arr.append(l2)

    l2_all = np.concatenate(l2_arr)
    n_snps = len(snp_arr)
    m_total = m_per_chrom * n_chrom

    # Trait 0: binary (with prevalences); trait 1: continuous.
    n1, n2 = 60_000.0, 50_000.0
    h2_1, h2_2, rg = 0.3, 0.4, 0.5
    infl1 = np.sqrt(1.0 + h2_1 * n1 * l2_all / m_total)
    infl2 = np.sqrt(1.0 + h2_2 * n2 * l2_all / m_total)
    g = rng.normal(0, 1, n_snps)
    e1 = rng.normal(0, 1, n_snps)
    e2 = rng.normal(0, 1, n_snps)
    z1 = infl1 * (np.sqrt(rg) * g + np.sqrt(1 - rg) * e1)
    z2 = infl2 * (np.sqrt(rg) * g + np.sqrt(1 - rg) * e2)

    # Random reference alleles so the flip logic is exercised.
    a1_t0 = rng.choice(["A", "G"], n_snps)
    a1_t1 = rng.choice(["A", "G"], n_snps)

    munged_paths = []
    for idx, (z, a1, n) in enumerate([(z1, a1_t0, n1), (z2, a1_t1, n2)]):
        df = pd.DataFrame(
            {
                "SNP": snp_arr,
                "A1": a1,
                "A2": np.where(a1 == "A", "G", "A"),
                "Z": z,
                "N": n,
            }
        )
        path = tmp_dir / f"trait{idx}.sumstats.gz"
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        munged_paths.append(path)

    return munged_paths, ld_dir, (0.3, float("nan")), (0.01, float("nan"))


@pytest.mark.skipif(
    not _have_rpy2_genomic_sem(), reason="rpy2/GenomicSEM not available"
)
def test_matches_genomic_sem_ldsc(tmp_path):
    """run_ldsc reproduces GenomicSEM::ldsc S, V, I, S_Stand, V_Stand."""
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    gsem = importr("GenomicSEM")

    n_chrom = 2
    munged_paths, ld_dir, samp_prev, pop_prev = _write_synthetic_inputs(
        tmp_path, n_chrom=n_chrom, snps_per_chrom=1500, seed=7
    )

    # ---- R ----
    r_result = gsem.ldsc(
        traits=ro.StrVector([str(p) for p in munged_paths]),
        sample_prev=ro.FloatVector(list(samp_prev)),
        population_prev=ro.FloatVector(list(pop_prev)),
        ld=str(ld_dir),
        wld=str(ld_dir),
        trait_names=ro.StrVector(["t0", "t1"]),
        n_blocks=200,
        chr=n_chrom,
        ldsc_log=str(tmp_path / "ldsc"),
        stand=True,
    )
    conv = ro.default_converter + pandas2ri.converter

    def r_mat(name: str) -> np.ndarray:
        with localconverter(conv):
            return np.asarray(ro.conversion.get_conversion().rpy2py(r_result.rx2(name)))

    r_S = r_mat("S")
    r_V = r_mat("V")
    r_I = r_mat("I")
    r_S_stand = r_mat("S_Stand")
    r_V_stand = r_mat("V_Stand")

    # ---- Python ----
    py = run_ldsc(
        munged_paths=munged_paths,
        ld_dir=ld_dir,
        sample_prev=list(samp_prev),
        population_prev=list(pop_prev),
        n_blocks=200,
        n_chrom=n_chrom,
        stand=True,
    )

    np.testing.assert_allclose(py.S, r_S, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose(py.I, r_I, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose(py.V, r_V, rtol=1e-4, atol=1e-10)
    assert py.S_Stand is not None and py.V_Stand is not None
    np.testing.assert_allclose(py.S_Stand, r_S_stand, rtol=1e-5, atol=1e-8)
    np.testing.assert_allclose(py.V_Stand, r_V_stand, rtol=1e-4, atol=1e-10)
