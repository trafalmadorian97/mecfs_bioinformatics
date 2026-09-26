"""
How each summary-statistic column changes when a variant's effect allele is swapped.

Columns are expected to use gwaslab's standard names (the "gwaslab" entry of formatbook.json), plus the regenie binary-trait columns and effective sample size used by datasets in this repo. Every non-allele column must be registered here, or declared through ExtraColumnRule, so that a new allele-dependent statistic can never be silently left unflipped.

All flipped values are computed from the original frame in a single with_columns call, so a confidence bound is never read after being overwritten.

FlipRule is what the type checker sees and the FLIP_* constants are what code uses. An import-time assertion keeps the two in step. A StrEnum cannot replace them, because polars converts Python Enum members into its own Enum dtype.
"""

from collections.abc import Mapping, Sequence
from typing import Final, Literal, get_args

import polars as pl
from attrs import frozen

from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHISQ_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_DOF_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_EFFECTIVE_SAMPLE_SIZE,
    GWASLAB_F_STATISTIC_COL,
    GWASLAB_HAZARD_RATIO_95L_COL,
    GWASLAB_HAZARD_RATIO_95U_COL,
    GWASLAB_HAZARD_RATIO_COL,
    GWASLAB_I2_COL,
    GWASLAB_INFO_SCORE_COL,
    GWASLAB_MAF_COL,
    GWASLAB_MLOG10P_COL,
    GWASLAB_N_CASE_COL,
    GWASLAB_N_CONTROL_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_ODDS_RATIO_95L_COL,
    GWASLAB_ODDS_RATIO_95U_COL,
    GWASLAB_ODDS_RATIO_COL,
    GWASLAB_P_COL,
    GWASLAB_P_HET_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
    GWASLAB_SNPID_COL,
    GWASLAB_SNPR2_COL,
    GWASLAB_STATUS_COL,
    GWASLAB_T_STATISTIC_COL,
    GWASLAB_Z_COL,
)
from mecfs_bio.constants.regenie_constants import (
    REGENIE_A1FREQ_CASES_COL,
    REGENIE_A1FREQ_CONTROLS_COL,
    REGENIE_EXTRA_COL,
    REGENIE_N_CASES_COL,
    REGENIE_N_CONTROLS_COL,
    REGENIE_TEST_COL,
)

FlipRule = Literal["negate", "complement", "invert", "invariant"]
FLIP_NEGATE: Final[FlipRule] = "negate"
FLIP_COMPLEMENT: Final[FlipRule] = "complement"
FLIP_INVERT: Final[FlipRule] = "invert"
FLIP_INVARIANT: Final[FlipRule] = "invariant"
assert set(get_args(FlipRule)) == {
    FLIP_NEGATE,
    FLIP_COMPLEMENT,
    FLIP_INVERT,
    FLIP_INVARIANT,
}


@frozen(slots=True)
class InvertedBoundPair:
    """Ratio confidence bounds: on a flip, lower becomes 1/upper and upper becomes 1/lower."""

    lower: str
    upper: str


@frozen(slots=True)
class ExtraColumnRule:
    column: str
    rule: FlipRule


COLUMN_FLIP_RULES: Mapping[str, FlipRule] = {
    GWASLAB_BETA_COL: FLIP_NEGATE,
    GWASLAB_Z_COL: FLIP_NEGATE,
    GWASLAB_T_STATISTIC_COL: FLIP_NEGATE,
    GWASLAB_EFFECT_ALLELE_FREQ_COL: FLIP_COMPLEMENT,
    REGENIE_A1FREQ_CASES_COL: FLIP_COMPLEMENT,
    REGENIE_A1FREQ_CONTROLS_COL: FLIP_COMPLEMENT,
    GWASLAB_ODDS_RATIO_COL: FLIP_INVERT,
    GWASLAB_HAZARD_RATIO_COL: FLIP_INVERT,
    GWASLAB_SNPID_COL: FLIP_INVARIANT,
    GWASLAB_RSID_COL: FLIP_INVARIANT,
    GWASLAB_CHROM_COL: FLIP_INVARIANT,
    GWASLAB_POS_COL: FLIP_INVARIANT,
    GWASLAB_SE_COL: FLIP_INVARIANT,
    GWASLAB_P_COL: FLIP_INVARIANT,
    GWASLAB_MLOG10P_COL: FLIP_INVARIANT,
    GWASLAB_CHISQ_COL: FLIP_INVARIANT,
    GWASLAB_F_STATISTIC_COL: FLIP_INVARIANT,
    GWASLAB_P_HET_COL: FLIP_INVARIANT,
    GWASLAB_I2_COL: FLIP_INVARIANT,
    GWASLAB_SNPR2_COL: FLIP_INVARIANT,
    GWASLAB_DOF_COL: FLIP_INVARIANT,
    GWASLAB_SAMPLE_SIZE_COLUMN: FLIP_INVARIANT,
    GWASLAB_N_CASE_COL: FLIP_INVARIANT,
    GWASLAB_N_CONTROL_COL: FLIP_INVARIANT,
    REGENIE_N_CASES_COL: FLIP_INVARIANT,
    REGENIE_N_CONTROLS_COL: FLIP_INVARIANT,
    GWASLAB_EFFECTIVE_SAMPLE_SIZE: FLIP_INVARIANT,
    GWASLAB_INFO_SCORE_COL: FLIP_INVARIANT,
    GWASLAB_MAF_COL: FLIP_INVARIANT,
    REGENIE_TEST_COL: FLIP_INVARIANT,
    REGENIE_EXTRA_COL: FLIP_INVARIANT,
}

BOUND_PAIRS: tuple[InvertedBoundPair, ...] = (
    InvertedBoundPair(
        lower=GWASLAB_ODDS_RATIO_95L_COL, upper=GWASLAB_ODDS_RATIO_95U_COL
    ),
    InvertedBoundPair(
        lower=GWASLAB_HAZARD_RATIO_95L_COL, upper=GWASLAB_HAZARD_RATIO_95U_COL
    ),
)

# Dropped from the output: STATUS describes gwaslab's processing, not this Task's.
DROPPED_COLUMNS: frozenset[str] = frozenset({GWASLAB_STATUS_COL})
_ALLELE_COLUMNS = frozenset({GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL})


def resolve_column_rules(
    columns: Sequence[str], extra: Sequence[ExtraColumnRule]
) -> Mapping[str, FlipRule]:
    """Flip rule for every present non-allele column; asserts that none is unregistered."""
    extra_rules = {rule.column: rule.rule for rule in extra}
    assert len(extra_rules) == len(extra), "extra_column_rules names a column twice"
    rules = {**COLUMN_FLIP_RULES, **extra_rules}
    bound_columns = {name for pair in BOUND_PAIRS for name in (pair.lower, pair.upper)}
    unregistered = [
        column
        for column in columns
        if column not in rules
        and column not in bound_columns
        and column not in _ALLELE_COLUMNS
        and column not in DROPPED_COLUMNS
    ]
    assert not unregistered, (
        f"columns without a flip rule: {unregistered}; register them in flip.py or pass "
        "ExtraColumnRule entries"
    )
    for pair in BOUND_PAIRS:
        assert (pair.lower in columns) == (pair.upper in columns), (
            f"confidence bounds {pair.lower} and {pair.upper} must be present together"
        )
    return {column: rules[column] for column in columns if column in rules}


def _flipped(column: str, rule: FlipRule) -> pl.Expr:
    value = pl.col(column)
    if rule == FLIP_NEGATE:
        return -value
    if rule == FLIP_COMPLEMENT:
        return 1 - value
    assert rule == FLIP_INVERT, f"rule {rule} has no flipped expression"
    return 1 / value


def flip_statistics(
    frame: pl.DataFrame, mask: pl.Expr, rules: Mapping[str, FlipRule]
) -> pl.DataFrame:
    """Apply every flip rule to the rows selected by mask; alleles are not touched."""
    updates = [
        pl.when(mask)
        .then(_flipped(column, rule))
        .otherwise(pl.col(column))
        .alias(column)
        for column, rule in rules.items()
        if rule != FLIP_INVARIANT and column in frame.columns
    ]
    for pair in BOUND_PAIRS:
        if pair.lower in frame.columns:
            updates.append(
                pl.when(mask)
                .then(1 / pl.col(pair.upper))
                .otherwise(pl.col(pair.lower))
                .alias(pair.lower)
            )
            updates.append(
                pl.when(mask)
                .then(1 / pl.col(pair.lower))
                .otherwise(pl.col(pair.upper))
                .alias(pair.upper)
            )
    return frame.with_columns(updates) if updates else frame
