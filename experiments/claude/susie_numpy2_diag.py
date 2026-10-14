"""
Diagnostic for the SUSIE adjustment discrepancy across the numpy 1.26 -> 2.5 upgrade.

Reproduces the numeric core of test_susie_r_finemap_task's failing case
(causal_variants = [], i.e. the null model) and fingerprints each intermediate so we
can localize where the two numpy versions diverge.

The raw PCG64 stream is guaranteed stable across numpy versions; the linear-algebra
transforms on top of it (multivariate_normal's covariance factorization, eigh in
make_psd_corr) are BLAS-dependent and are the suspects.

Run on each branch/env:
    pixi r python experiments/claude/susie_numpy2_diag.py 2>&1 | tee experiments/claude/logs/susie_numpy2_diag_<label>.log
"""

import hashlib

import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

# Import the REAL production helper so the eigh-based PSD step matches the task exactly.
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import make_psd_corr


def fp(name: str, a) -> None:
    arr = np.ascontiguousarray(np.asarray(a, dtype=np.float64))
    print(
        f"[DIAG] {name:22s} shape={str(arr.shape):12s} "
        f"sum={arr.sum():+.12e} sq={np.square(arr).sum():+.12e} "
        f"min={arr.min():+.6e} max={arr.max():+.6e} "
        f"md5={hashlib.md5(arr.tobytes()).hexdigest()}",
        flush=True,
    )


def main() -> None:
    print(f"numpy version: {np.__version__}", flush=True)
    print("numpy BLAS config:", flush=True)
    np.show_config()

    n = 2500
    m = 100
    susie_package = importr("susieR")

    # --- Prove the raw RNG stream is version-independent (should match across envs). ---
    fp("raw_normals", np.random.default_rng(40).standard_normal(2000))

    # --- Reproduce the fixture (causal_variants = [], the failing null case). ---
    generator = np.random.default_rng(40)
    true_effects = np.zeros(m)  # causal_variants = [] -> all zero
    q = np.linspace(1, 0, num=m)
    lamb = 0.2
    covar = (1 - lamb) * q.reshape(-1, 1) * q.reshape(1, -1) + lamb * np.eye(m)
    fp("covar", covar)

    genotypes = generator.multivariate_normal(np.zeros(m), cov=covar, size=n)
    fp("genotypes", genotypes)

    phenotypes = genotypes @ true_effects.reshape(-1) + generator.normal(
        loc=0, scale=0.1, size=n
    )
    fp("phenotypes", phenotypes)

    corr = np.corrcoef(genotypes.transpose())
    fp("corrcoef", corr)

    # The task's LD input is corrcoef / 2 (see the fixture's partial_ld); permutation
    # (the flip) does not change estimate_s_rss, so we omit it.
    partial_ld = corr / 2
    fp("partial_ld", partial_ld)

    conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter
    with localconverter(conv):
        genotypes_r = ro.conversion.get_conversion().py2rpy(genotypes)
        phenotypes_r = ro.conversion.get_conversion().py2rpy(phenotypes)
    reg = susie_package.univariate_regression(genotypes_r, phenotypes_r)
    with localconverter(conv):
        beta_hat = np.asarray(reg.rx2("betahat"))
        se_beta_hat = np.asarray(reg.rx2("sebetahat"))
    fp("beta_hat", beta_hat)
    fp("se_beta_hat", se_beta_hat)

    zscores = beta_hat / se_beta_hat
    fp("zscores", zscores)

    # Task path: make_psd_corr (eigh) then estimate_s_rss.
    eigs = np.linalg.eigvalsh(partial_ld)
    print(
        f"[DIAG] eigvals of partial_ld: min={eigs.min():+.6e} "
        f"n_below_tol(1e-4)={(eigs < 1e-4).sum()}",
        flush=True,
    )
    ld_psd = make_psd_corr(partial_ld)
    fp("ld_after_make_psd", ld_psd)

    with localconverter(conv):
        zscores_r = ro.conversion.get_conversion().py2rpy(zscores)
        ld_r = ro.conversion.get_conversion().py2rpy(ld_psd)
        adjustment = susie_package.estimate_s_rss(zscores_r, ld_r, n=n)
    adj = float(np.asarray(adjustment).item())
    print(f"[DIAG] >>> ADJUSTMENT = {adj:.12e}  (test bound: <= 0.01)", flush=True)


if __name__ == "__main__":
    main()
