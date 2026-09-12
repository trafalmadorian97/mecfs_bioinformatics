"""
Render a realistic synthetic LD-score-regression diagnostic plot to eyeball the design.

Mimics a DecodeME-like case: an intercept slightly below 1 (the confounding signal the plot is
meant to surface) plus a real heritability slope. Not a test -- just a visual check.
Run: pixi r python experiments/claude/ldsc_diagnostic/preview.py
"""

from pathlib import Path

import numpy as np

from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic import bin_by_ld_score
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic_plot_task import (
    LdscDiagnosticPlotConfig,
    LdscFit,
    build_diagnostic_figure,
)

OUT = Path(__file__).parent
rng = np.random.default_rng(0)

n = 200_000.0
m = 1_000_000.0
true_h2 = 0.2
true_intercept = 0.93  # below 1, like the DecodeME diagnostic
slope = n * true_h2 / m

n_snps = 200_000
ld = rng.gamma(shape=2.0, scale=40.0, size=n_snps)  # mean ~80, right-skewed like real LD scores
expected = true_intercept + slope * ld
chi2 = expected * rng.chisquare(df=1, size=n_snps)  # E[chi^2_1df] = 1

# In production the fit comes from gwaslab; here we hand it the true values to eyeball the figure.
fit = LdscFit(intercept=true_intercept, h2_obs=true_h2)
bins = bin_by_ld_score(chi2, ld, n_bins=25)

fig = build_diagnostic_figure(
    bins=bins,
    fit=fit,
    n=n,
    m=m,
    config=LdscDiagnosticPlotConfig(
        n_bins=25, show_error_bars=True, title="LDSC diagnostic (synthetic DecodeME-like)"
    ),
)
fig.write_html(str(OUT / "preview.html"), include_plotlyjs="cdn")
try:
    fig.write_image(str(OUT / "preview.png"), width=900, height=600, scale=2)
    print(f"wrote {OUT / 'preview.png'}")
except Exception as exc:  # kaleido may be absent
    print(f"PNG export unavailable ({exc}); wrote HTML only")
