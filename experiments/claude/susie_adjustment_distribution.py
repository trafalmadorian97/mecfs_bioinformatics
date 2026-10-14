"""
Monte-Carlo the SUSIE LD-adjustment (estimate_s_rss) over many random draws of the
test's null-model data-generating process, to calibrate a distribution-level upper bound
for the test that does NOT depend on a specific realization.

Faithful to the real task: the LD fed to estimate_s_rss is the FULL correlation matrix
(the fixture stores corrcoef/2 and the task symmetrizes via M + M.T -> corrcoef, diag 1),
then passed through the production make_psd_corr.

Run:
    pixi r python experiments/claude/susie_adjustment_distribution.py 2>&1 | tee experiments/claude/logs/susie_adjustment_distribution.log
"""

import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import make_psd_corr

M = 100
N = 2500
LAMB = 0.2
_susie = importr("susieR")
_conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter


def adjustment_for_seed(seed: int) -> float:
    gen = np.random.default_rng(seed)
    q = np.linspace(1, 0, num=M)
    covar = (1 - LAMB) * q.reshape(-1, 1) * q.reshape(1, -1) + LAMB * np.eye(M)
    genotypes = gen.multivariate_normal(np.zeros(M), cov=covar, size=N)
    # Null model: true_effects = 0, so phenotype is pure noise.
    phenotypes = gen.normal(loc=0, scale=0.1, size=N)

    corr = np.corrcoef(genotypes.T)  # diag 1, matches task after M + M.T
    with localconverter(_conv):
        genotypes_r = ro.conversion.get_conversion().py2rpy(genotypes)
        phenotypes_r = ro.conversion.get_conversion().py2rpy(phenotypes)
    reg = _susie.univariate_regression(genotypes_r, phenotypes_r)
    with localconverter(_conv):
        beta = np.asarray(reg.rx2("betahat"))
        se = np.asarray(reg.rx2("sebetahat"))
    z = beta / se
    ld = make_psd_corr(corr)
    with localconverter(_conv):
        adj = _susie.estimate_s_rss(
            ro.conversion.get_conversion().py2rpy(z),
            ro.conversion.get_conversion().py2rpy(ld),
            n=N,
        )
    return float(np.asarray(adj).item())


def main() -> None:
    print(f"numpy {np.__version__}", flush=True)
    print("validate: seed 40 adjustment =", f"{adjustment_for_seed(40):.6f}", flush=True)

    seeds = range(1000, 1000 + 300)
    vals = []
    for s in seeds:
        vals.append(adjustment_for_seed(s))
        if len(vals) % 25 == 0:
            print(f"  ...{len(vals)} draws done", flush=True)
    a = np.array(vals)
    qs = [50, 90, 95, 99, 99.9, 100]
    print(f"\nn_draws = {a.size}")
    print(f"mean    = {a.mean():.6f}")
    print(f"std     = {a.std():.6f}")
    for qv in qs:
        print(f"p{qv:<5} = {np.percentile(a, qv):.6f}")
    print(f"min/max = {a.min():.6f} / {a.max():.6f}")


if __name__ == "__main__":
    main()
