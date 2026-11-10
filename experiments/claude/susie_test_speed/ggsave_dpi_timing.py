"""
Time susieR::kriging_rss and ggsave of its diagnostic plot at several dpi values,
for an m=100 locus like the SuSiE unit-test fixture.

Usage: pixi r python experiments/claude/susie_test_speed/ggsave_dpi_timing.py
"""

import tempfile
import time
from pathlib import Path

import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter

from mecfs_bio.util.r_packages.cached_importr import cached_importr

susie = cached_importr("susieR")
ggplot2 = cached_importr("ggplot2")
m = 100
q = np.linspace(1, 0, num=m)
ld = 0.8 * np.outer(q, q) + 0.2 * np.eye(m)
z = np.random.default_rng(0).normal(size=m)
with localconverter(ro.default_converter + numpy2ri.converter):
    z_r = ro.conversion.get_conversion().py2rpy(z)
    ld_r = ro.conversion.get_conversion().py2rpy(ld)
for rep in range(3):
    t = time.perf_counter()
    res = susie.kriging_rss(z_r, ld_r, n=2500)
    print(f"kriging_rss {time.perf_counter() - t:.3f}s")
plot = res.rx2("plot")
with tempfile.TemporaryDirectory() as d:
    for dpi in [300, 300, 100, 72]:
        t = time.perf_counter()
        ggplot2.ggsave(filename=str(Path(d) / "p.png"), plot=plot, width=7, height=7, dpi=dpi, device="png")
        print(f"ggsave dpi={dpi} {time.perf_counter() - t:.3f}s")
