"""
Script to compute all-pairs genetic correlation for the initial study traits with the fast
CT-LDSC generator, and plot the result.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_fast_initial_asset_generator import (
    CT_LDSC_FAST_INITIAL_ASSET_GENERATOR,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_fast_plot import (
    CT_LDSC_FAST_INITIAL_PLOT,
)


def fast_ct_ldsc_analysis() -> None:
    """
    Compute fast all-pairs CT-LDSC genetic correlation and heritability for the initial traits.
    """
    DEFAULT_RUNNER.run(
        CT_LDSC_FAST_INITIAL_ASSET_GENERATOR.terminal_tasks()
        + [CT_LDSC_FAST_INITIAL_PLOT],
        incremental_save=True,
    )


if __name__ == "__main__":
    fast_ct_ldsc_analysis()
