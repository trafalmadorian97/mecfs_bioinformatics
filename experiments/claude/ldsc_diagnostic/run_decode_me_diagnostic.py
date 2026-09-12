"""
Drive the DecodeME LD-score-regression diagnostic plot task and build it.

Run: pixi r python experiments/claude/ldsc_diagnostic/run_decode_me_diagnostic.py
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_ldsc_diagnostic_plot import (
    DECODE_ME_GWAS_1_LDSC_DIAGNOSTIC_PLOT,
)


def go():
    DEFAULT_RUNNER.run(targets=[DECODE_ME_GWAS_1_LDSC_DIAGNOSTIC_PLOT], must_rebuild_transitive=[DECODE_ME_GWAS_1_LDSC_DIAGNOSTIC_PLOT])


if __name__ == "__main__":
    go()
