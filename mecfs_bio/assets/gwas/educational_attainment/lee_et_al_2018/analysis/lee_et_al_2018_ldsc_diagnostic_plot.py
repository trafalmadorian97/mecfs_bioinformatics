"""
LD-score-regression diagnostic plot for the Lee et al. 2018 educational-attainment heritability run.

Bins the variants by LD score and shows each bin's mean chi-square against the fitted line, so model
misfit and confounding in the heritability estimate can be inspected by eye. The fit is re-run with
gwaslab -- the same estimator lee_et_al_h2_by_ldsc uses -- so the line and annotation match that run.
"""

from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.analysis.lee_et_al_h2_by_ldsc import (
    LEE_ET_AL_2018_H2_BY_LDSC,
)
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic_plot_task import (
    LdscDiagnosticPlotConfig,
    LdscDiagnosticPlotTask,
)

LEE_ET_AL_2018_LDSC_DIAGNOSTIC_PLOT = LdscDiagnosticPlotTask.create(
    asset_id="lee_et_al_2018_ldsc_diagnostic_plot",
    ldsc_task=LEE_ET_AL_2018_H2_BY_LDSC,
    config=LdscDiagnosticPlotConfig(
        n_bins=25,
        show_error_bars=True,
        title="Lee et al. 2018 educational attainment — LDSC diagnostic",
    ),
)
