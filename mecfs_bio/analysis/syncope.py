"""
Analyze syncope GWAS GWAS
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.analysis.syncope_s_ldsc import (
    SYNCOPE_S_LDSC_TASKS,
)
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_magma_task_generator import (
    AEGISDOTTIR_COMBINED_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.raw.raw_syncope_data import (
    AEGISDOTTIR_ET_AL_RAW_SYNCOPE_DATA,
)


def run_syncope_analysis():
    """
    Analyze GWAS of syncope from aegisdottir et al.


    """
    DEFAULT_RUNNER.run(
        [
            AEGISDOTTIR_ET_AL_RAW_SYNCOPE_DATA,
        ]
        + SYNCOPE_S_LDSC_TASKS.get_terminal_tasks(),
        incremental_save=True,
    )


if __name__ == "__main__":
    run_syncope_analysis()
