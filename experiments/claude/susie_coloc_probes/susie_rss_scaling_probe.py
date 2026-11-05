"""
Timing probe: cost of loading a Broad UKBB LD matrix (.npz, lower triangle with halved diagonal)
and of calling R susie_rss through rpy2 at coloc-sized windows.

Usage: pixi r python susie_rss_scaling_probe.py <path to Broad .npz>

Takes the first p variants of the matrix for several p, simulates z-scores with two causal
variants, z = sqrt(N) R b + R^(1/2) e with e standard normal and negative eigenvalues of R
clipped to zero, and times the R-to-numpy
conversion plus susie_rss with L=10.
"""

import sys
import time
from pathlib import Path

import numpy as np
import rpy2.robjects as ro
import scipy.sparse
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

N = 50_000


def main(npz_path: Path) -> None:
    t0 = time.perf_counter()
    lower = scipy.sparse.load_npz(npz_path).tocsr()
    print(f"load_npz: {time.perf_counter() - t0:.1f} s, shape={lower.shape}, nnz={lower.nnz}")
    susie = importr("susieR")
    rng = np.random.default_rng(0)
    for p in (1000, 2000, 4000, 6000):
        t0 = time.perf_counter()
        block = lower[:p, :p].toarray()
        ld = block + block.T
        t_dense = time.perf_counter() - t0
        b = np.zeros(p)
        b[[p // 3, 2 * p // 3]] = [0.04, -0.03]
        evals, evecs = np.linalg.eigh(ld)
        noise = evecs @ (np.sqrt(np.clip(evals, 0.0, None)) * rng.standard_normal(p))
        z = ld @ b * np.sqrt(N) + noise
        t0 = time.perf_counter()
        with localconverter(ro.default_converter + numpy2ri.converter):
            ld_r = ro.conversion.get_conversion().py2rpy(ld)
            z_r = ro.conversion.get_conversion().py2rpy(z)
        t_conv = time.perf_counter() - t0
        t0 = time.perf_counter()
        fit = susie.susie_rss(z=z_r, R=ld_r, n=N, L=10)
        t_fit = time.perf_counter() - t0
        niter = int(np.asarray(fit.rx2("niter"))[0])
        print(f"p={p}: densify {t_dense:.2f} s, py2rpy {t_conv:.2f} s, susie_rss {t_fit:.2f} s ({niter} iters)")
        del ld_r, fit
        ro.r("gc()")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
