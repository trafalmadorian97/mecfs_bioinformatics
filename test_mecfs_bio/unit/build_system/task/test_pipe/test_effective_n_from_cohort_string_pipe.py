import narwhals
import pandas as pd

from mecfs_bio.build_system.task.pipes.effective_n_from_cohort_string_pipe import (
    CohortCaseControl,
    EffectiveNFromCohortStringPipe,
)


def _round_eff(n_cases: int, n_controls: int) -> int:
    return round(4.0 * n_cases * n_controls / (n_cases + n_controls))


def test_effective_n_sums_contributing_cohorts():
    # Two cohorts in fixed order; position 0 contributes when char 0 != "?", etc.
    cohorts = [
        CohortCaseControl(name="A", n_cases=1000, n_controls=4000),
        CohortCaseControl(name="B", n_cases=500, n_controls=500),
    ]
    eff_a = _round_eff(1000, 4000)
    eff_b = _round_eff(500, 500)

    data = pd.DataFrame({"Cohorts": ["+?", "?-", "+-", "??"]})
    pipe = EffectiveNFromCohortStringPipe(
        cohorts=cohorts, cohort_string_col="Cohorts", out_col="N"
    )
    result = pipe.process(narwhals.from_native(data).lazy()).collect().to_pandas()

    assert result["N"].tolist() == [eff_a, eff_b, eff_a + eff_b, 0]


def test_effective_n_of_cohort_case_control():
    # Balanced 1:1 cohort: effective N equals the total sample size.
    assert CohortCaseControl(name="X", n_cases=100, n_controls=100).effective_n == 200.0
