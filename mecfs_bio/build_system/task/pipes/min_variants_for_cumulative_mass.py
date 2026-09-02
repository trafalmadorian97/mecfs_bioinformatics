"""Pipe: within each group, keep the smallest set of highest-value rows whose
cumulative value first reaches a threshold."""

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe

_NEG_COL = "__neg_value"
_PREV_CUM_COL = "__prev_cumulative"


@frozen
class MinVariantsForCumulativeMass(DataProcessingPipe):
    """Within each group_col, keep the minimal prefix of highest-value_col rows
    whose cumulative value_col first reaches threshold (the crossing row is
    included), dropping the low-value tail. Applied to a SUSIE credible-set table
    (group_col = credible set, value_col = PIP) this is the per-credible-set
    threshold% credible set, unioned across sets.

    A single row with value >= threshold yields just that row; a group whose
    total never reaches threshold keeps all its rows.

    Implemented for a lazy frame: cumulative sum needs an explicit order, so the
    rows are ordered by descending value_col (via a negated helper column, since
    the order is ascending), and a row is kept when the cumulative value of the
    strictly-higher rows before it is still below threshold."""

    group_col: str
    value_col: str
    threshold: float = 0.5

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return (
            x.with_columns((-narwhals.col(self.value_col)).alias(_NEG_COL))
            .with_columns(
                (
                    narwhals.col(self.value_col)
                    .cum_sum()
                    .over(self.group_col, order_by=_NEG_COL)
                    - narwhals.col(self.value_col)
                ).alias(_PREV_CUM_COL)
            )
            .filter(narwhals.col(_PREV_CUM_COL) < self.threshold)
            .drop(_NEG_COL, _PREV_CUM_COL)
        )
