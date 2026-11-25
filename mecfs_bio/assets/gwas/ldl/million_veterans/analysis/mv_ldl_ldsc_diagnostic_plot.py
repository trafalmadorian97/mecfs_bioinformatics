"""
LD-score-regression diagnostic plot for the Million Veterans Program LDL heritability run.

Bins the LDL variants by LD score and shows each bin's mean chi-square against the fitted line, so
model misfit and confounding in the heritability estimate can be inspected by eye. The fit is
re-run with gwaslab -- the same estimator mv_ldl_heritability_task uses -- so the line and
annotation match that run.
"""

from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.mv_ldl_heritability_task import (
    MV_LDL_HERITABILITY_TASK,
)
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic_plot_task import (
    LdscDiagnosticPlotConfig,
    LdscDiagnosticPlotTask,
)

MV_LDL_LDSC_DIAGNOSTIC_PLOT = LdscDiagnosticPlotTask.create(
    asset_id="million_veterans_eur_ldl_ldsc_diagnostic_plot",
    ldsc_task=MV_LDL_HERITABILITY_TASK,
    config=LdscDiagnosticPlotConfig(
        n_bins=25,
        show_error_bars=True,
        title="Million Veterans LDL — LDSC diagnostic",
    ),
)
