"""
Script to perform standard analysis on the FinnGen release 13 GWAS of migraine and its subtypes.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.migraine.finngen_r13.analysis.finngen_r13_migraine_standard_analysis import (
    FINNGEN_R13_MIGRAINE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_no_aura.analysis.finngen_r13_migraine_no_aura_standard_analysis import (
    FINNGEN_R13_MIGRAINE_NO_AURA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_triptan.analysis.finngen_r13_migraine_triptan_standard_analysis import (
    FINNGEN_R13_MIGRAINE_TRIPTAN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.migraine.finngen_r13_migraine_with_aura.analysis.finngen_r13_migraine_with_aura_standard_analysis import (
    FINNGEN_R13_MIGRAINE_WITH_AURA_STANDARD_ANALYSIS,
)


def run_finngen_r13_migraine_analysis():
    """
    Function to perform standard analysis on the FinnGen R13 migraine endpoints:
    G6_MIGRAINE, MIGRAINE_TRIPTAN, G6_MIGRAINE_NO_AURA, and G6_MIGRAINE_WITH_AURA.

    Includes S-LDSC, LDSC heritability, GTEx MAGMA, HBA MAGMA, H-MAGMA, gene set analysis,
    lead variant extraction, and Manhattan plots.
    """
    DEFAULT_RUNNER.run(
        FINNGEN_R13_MIGRAINE_STANDARD_ANALYSIS.get_terminal_tasks()
        + FINNGEN_R13_MIGRAINE_TRIPTAN_STANDARD_ANALYSIS.get_terminal_tasks()
        + FINNGEN_R13_MIGRAINE_NO_AURA_STANDARD_ANALYSIS.get_terminal_tasks()
        + FINNGEN_R13_MIGRAINE_WITH_AURA_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_finngen_r13_migraine_analysis()
