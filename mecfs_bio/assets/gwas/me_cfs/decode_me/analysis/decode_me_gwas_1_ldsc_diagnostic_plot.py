"""
LD-score-regression diagnostic plot for the main DecodeME LDSC run.

Bins the DecodeME variants by LD score and shows each bin's mean chi-square against the fitted
line, so the intercept-below-one that the heritability run reports can be inspected directly (does
the low-LD end sit below one? do the bins fall on the line?). The fit is re-run with gwaslab -- the
same estimator decode_me_gwas_1_ldsc uses -- so the line and annotation match that run.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_ldsc import (
    DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC,
)
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic_plot_task import (
    LdscDiagnosticPlotConfig,
    LdscDiagnosticPlotTask,
)

DECODE_ME_GWAS_1_LDSC_DIAGNOSTIC_PLOT = LdscDiagnosticPlotTask.create(
    asset_id="decode_me_gwas_1_ldsc_diagnostic_plot",
    ldsc_task=DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC,
    config=LdscDiagnosticPlotConfig(
        n_bins=25,
        show_error_bars=True,
        title="DecodeME GWAS 1 — LDSC diagnostic",
    ),
)
