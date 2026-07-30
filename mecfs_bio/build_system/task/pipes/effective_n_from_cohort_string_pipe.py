"""
Pipe that derives a per-variant effective sample size for a case/control
meta-analysis whose summary statistics encode, per variant, which cohorts
contributed via a fixed-order string of per-cohort characters.

The motivating example is the deCODE Cohorts column, where each character is +, -
or ?, one position per cohort in a fixed order, and ? means the cohort did not
contribute to that variant. Different variants are therefore backed by different
subsets of cohorts and so have different sample sizes, even though the file has no
sample-size column.

For a case/control meta-analysis the sample size that downstream tools (LDSC,
MAGMA) want is the sum over contributing cohorts of each cohort's effective sample
size, 4 / (1/n_cases + 1/n_controls). Summing per-cohort effective N, rather than
assigning the whole study's total N to every variant, correctly down-weights
variants that only a subset of cohorts contributed to and accounts for
case/control imbalance within each cohort.
"""

import functools
import operator

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN


@frozen
class CohortCaseControl:
    """Case and control counts for a single cohort of a case/control meta-analysis."""

    name: str
    n_cases: int
    n_controls: int

    def __attrs_post_init__(self) -> None:
        assert self.n_cases > 0, f"cohort {self.name} has non-positive case count"
        assert self.n_controls > 0, f"cohort {self.name} has non-positive control count"

    @property
    def effective_n(self) -> float:
        """Effective sample size 4 / (1/n_cases + 1/n_controls)."""
        return 4.0 * self.n_cases * self.n_controls / (self.n_cases + self.n_controls)


@frozen
class EffectiveNFromCohortStringPipe(DataProcessingPipe):
    """Add a per-variant effective sample size column derived from a per-variant
    cohort-membership string.

    cohorts is ordered to match the character positions of cohort_string_col:
    element i contributes its effective N to a variant exactly when character i of
    that variant's string is not missing_char.
    """

    cohorts: list[CohortCaseControl]
    cohort_string_col: str
    out_col: str = GWASLAB_SAMPLE_SIZE_COLUMN
    missing_char: str = "?"

    def __attrs_post_init__(self) -> None:
        assert len(self.missing_char) == 1, "missing_char must be a single character"
        assert self.cohorts, "at least one cohort is required"

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        terms = [
            narwhals.when(
                narwhals.col(self.cohort_string_col).str.slice(position, 1)
                != self.missing_char
            )
            .then(narwhals.lit(cohort.effective_n))
            .otherwise(narwhals.lit(0.0))
            for position, cohort in enumerate(self.cohorts)
        ]
        effective_n = functools.reduce(operator.add, terms)
        return x.with_columns(
            effective_n.round(0).cast(narwhals.Int64).alias(self.out_col)
        )
