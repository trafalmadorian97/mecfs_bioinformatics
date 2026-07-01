"""
Script to generate a heatmap showing estimated genetic correlation between a diverse collection of traits.
Uses cross-trait linkage disequilibrium score regression to estimate correlation.

See:
Bulik-Sullivan, Brendan, et al. "An atlas of genetic correlations across human diseases and traits." Nature genetics 47.11 (2015): 1236-1241.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_initial_asset_generator import (
    CT_LDSC_INITIAL_ASSET_GENERATOR,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_plot import (
    CT_LDSC_INITIAL_PLOT,
)
from mecfs_bio.build_system.scheduler.topological_scheduler import (
    TopologicalSchedulerSettings,
)


def initial_ct_ldsc_analysis():
    DEFAULT_RUNNER.run(
        [CT_LDSC_INITIAL_PLOT],
        incremental_save=True,
        must_rebuild_transitive=[CT_LDSC_INITIAL_ASSET_GENERATOR.aggregation_task],
        settings=TopologicalSchedulerSettings(print_progress=True),
    )


if __name__ == "__main__":
    initial_ct_ldsc_analysis()
