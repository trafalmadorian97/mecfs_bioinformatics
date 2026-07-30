"""
Script to use CT-LDSC to compute genetic correlation between various CFS-related GWAS
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.multistudy.analysis.genetic_correlation.ct_ldsc.ct_ldsc_mecfs_studies import (
    CFS_CT_LDSC_ASSET_GENERATOR,
)
from mecfs_bio.assets.gwas.me_cfs.multistudy.analysis.genetic_correlation.ct_ldsc.ct_ldsc_mecfs_studies_plot import (
    CT_LDSC_CFS_CORR_PLOT,
)


def cfs_genetic_corr():
    """
    Compute genetic correlation between CFS-related GWAS
    Plot the result.
    """
    DEFAULT_RUNNER.run(
        (CFS_CT_LDSC_ASSET_GENERATOR.terminal_tasks() + [CT_LDSC_CFS_CORR_PLOT]),
        incremental_save=True,
        must_rebuild_transitive=[CT_LDSC_CFS_CORR_PLOT],
    )


if __name__ == "__main__":
    cfs_genetic_corr()
