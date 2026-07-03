"""
Loader for the per-cohort case/control counts of the deCODE rheumatoid arthritis
meta-analyses (Saevarsdottir et al 2022, Ann Rheum Dis 81:1085-1095, Table 1).

The counts live in decode_cohort_sample_sizes.csv beside this module; see the
accompanying decode_cohort_sample_sizes.md for provenance. They are used to turn
the per-variant deCODE Cohorts membership string into a per-variant effective
sample size via EffectiveNFromCohortStringPipe.
"""

from pathlib import Path
from typing import Literal

import polars as pl

from mecfs_bio.build_system.task.pipes.effective_n_from_cohort_string_pipe import (
    CohortCaseControl,
)

Serostatus = Literal["seropositive", "seronegative"]

_COHORT_TABLE_PATH = Path(
    "mecfs_bio/assets/gwas/rheumtoid_arthritis/util/decode_cohort_sample_sizes.csv"
)


def load_decode_ra_cohorts(serostatus: Serostatus) -> list[CohortCaseControl]:
    """Cohorts for the given serostatus, ordered to match the character positions
    of the deCODE Cohorts string (Iceland, Norway, Sweden, Denmark, UK, Finland)."""
    table = pl.read_csv(_COHORT_TABLE_PATH).sort("position")
    return [
        CohortCaseControl(
            name=row["cohort"],
            n_cases=row[f"{serostatus}_cases"],
            n_controls=row[f"{serostatus}_controls"],
        )
        for row in table.iter_rows(named=True)
    ]
