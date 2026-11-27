"""
Driver to colocalize traits with the deCODE SomaScan SMP assays: raw results and display table
per run. The per-protein deCODE database must be built first
(decode_somascan_database_build.py), since each run reads every protein's stored seek index.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_decode_somascan_coloc_exploratory import (
    DECODE_ME_DECODE_SOMASCAN_SMP_COLOC_LEAD_P_1E_5,
    DECODE_ME_DECODE_SOMASCAN_SMP_COLOC_LEAD_P_1E_5_PLEIOTROPY,
)


def run_decode_me_exploratory() -> None:
    DEFAULT_RUNNER.run(
        [
            *DECODE_ME_DECODE_SOMASCAN_SMP_COLOC_LEAD_P_1E_5.terminal_tasks(),
            DECODE_ME_DECODE_SOMASCAN_SMP_COLOC_LEAD_P_1E_5_PLEIOTROPY,
        ],
        incremental_save=True,
        must_rebuild_transitive=[],
    )


def main() -> None:
    run_decode_me_exploratory()


if __name__ == "__main__":
    main()
