"""
Why is the estimate_s_rss adjustment sometimes large when z and R come from the SAME data?

Hypothesis: the null-mle for s is dominated by the small-eigenvalue directions of R (the
LD matrix is ill-conditioned because the SNPs share a strong common factor q, so many
eigenvalues are tiny).  In eigenbasis u = V^T z with eigenvalues d, under correct
specification u_k ~ N(0, d_k), so the leverage u_k^2 / d_k ~ 1 in expectation but has
variance ~ 1/d_k.  When a random z projects heavily onto a tiny-d_k direction, u_k^2/d_k
blows up and the null-mle inflates s to add variance there.  So s is essentially a
high-variance diagnostic of z's projection onto R's smallest-eigenvalue directions, NOT a
sign that z and R systematically disagree (the median s is 0).

Run:
    pixi r python experiments/claude/susie_adjustment_mechanism.py 2>&1 | tee experiments/claude/logs/susie_adjustment_mechanism.log
"""

import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from scipy.stats import spearmanr

M, N, LAMB = 100, 2500, 0.2
_susie = importr("susieR")
_conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter
_Q = np.linspace(1, 0, num=M)
_COVAR = (1 - LAMB) * _Q.reshape(-1, 1) * _Q.reshape(1, -1) + LAMB * np.eye(M)


def draw(seed: int):
    gen = np.random.default_rng(seed)
    genotypes = gen.multivariate_normal(np.zeros(M), cov=_COVAR, size=N)
    phenotypes = gen.normal(loc=0, scale=0.1, size=N)  # null: y independent of genotypes
    corr = np.corrcoef(genotypes.T)
    with localconverter(_conv):
        g_r = ro.conversion.get_conversion().py2rpy(genotypes)
        p_r = ro.conversion.get_conversion().py2rpy(phenotypes)
    reg = _susie.univariate_regression(g_r, p_r)
    with localconverter(_conv):
        beta = np.asarray(reg.rx2("betahat"))
        se = np.asarray(reg.rx2("sebetahat"))
    z = beta / se
    with localconverter(_conv):
        adj = _susie.estimate_s_rss(
            ro.conversion.get_conversion().py2rpy(z),
            ro.conversion.get_conversion().py2rpy(corr),
            n=N,
        )
    s = float(np.asarray(adj).item())
    return z, corr, s


def leverage_stats(z, corr):
    d, V = np.linalg.eigh(corr)  # ascending eigenvalues
    u = V.T @ z
    lev = u**2 / d  # per-direction chi-square-like leverage
    return d, lev, u


def main() -> None:
    print(f"numpy {np.__version__}\n", flush=True)

    # 1. Show R's eigenspectrum for a representative draw.
    z0, corr0, s0 = draw(40)
    d0, lev0, u0 = leverage_stats(z0, corr0)
    print("Representative draw (seed 40, s = %.4f):" % s0)
    print("  R eigenvalues: min=%.3e  median=%.3e  max=%.3e  cond=%.3e"
          % (d0.min(), np.median(d0), d0.max(), d0.max() / d0.min()))
    order = np.argsort(lev0)[::-1]
    print("  top-5 leverage directions (u_k^2/d_k):")
    for k in order[:5]:
        print(f"    d_k={d0[k]:.3e}  u_k^2={u0[k] ** 2:.3e}  leverage={lev0[k]:.2f}")

    # 2. Across many draws, correlate s with the ill-conditioning / leverage signals.
    seeds = range(2000, 2000 + 200)
    S, MINEIG, MAXLEV, LEV_SMALL = [], [], [], []
    for sd in seeds:
        z, corr, s = draw(sd)
        d, lev, _u = leverage_stats(z, corr)
        S.append(s)
        MINEIG.append(d.min())
        MAXLEV.append(lev.max())
        # summed leverage in the 10 smallest-eigenvalue directions
        LEV_SMALL.append(lev[np.argsort(d)[:10]].sum())
    S, MINEIG, MAXLEV, LEV_SMALL = map(np.array, (S, MINEIG, MAXLEV, LEV_SMALL))

    pos = S > 0
    print(f"\nOver {S.size} draws:  frac(s==0) = {(~pos).mean():.2f}  "
          f"mean_s = {S.mean():.4f}  max_s = {S.max():.4f}")
    print("Spearman correlations with s:")
    print(f"  max leverage max_k u_k^2/d_k         : rho = {spearmanr(S, MAXLEV).statistic:+.3f}")
    print(f"  summed leverage in 10 smallest dirs   : rho = {spearmanr(S, LEV_SMALL).statistic:+.3f}")
    print(f"  min eigenvalue of R (smaller->bigger s?): rho = {spearmanr(S, MINEIG).statistic:+.3f}")
    print("\n(A strong positive rho with leverage confirms: large s is driven by z projecting"
          "\n onto R's smallest-eigenvalue directions, not by a systematic z/R mismatch.)")


if __name__ == "__main__":
    main()
