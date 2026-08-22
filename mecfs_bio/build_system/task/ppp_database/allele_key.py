"""
Order-agnostic allele key, shared across the PPP-database tasks for allele-aware
variant matching (so a swapped effect/non-effect orientation still matches).
"""

import narwhals as nw
import polars as pl


def unordered_allele_key(a: str, b: str) -> pl.Expr:
    """Sort the two allele columns and join them, so {A, B} == {B, A}."""
    return (
        pl.when(pl.col(a) <= pl.col(b))
        .then(pl.col(a) + pl.lit("_") + pl.col(b))
        .otherwise(pl.col(b) + pl.lit("_") + pl.col(a))
    )


def unordered_allele_key_narwhals(a: str, b: str) -> nw.Expr:
    """Narwhals sibling of unordered_allele_key, for building the key while staying
    on a narwhals frame (so the whole pipeline defers to a single collect).

    Kept identical to the polars version above; the two are exercised against each
    other in the tests so they cannot silently diverge.
    """
    return (
        nw.when(nw.col(a) <= nw.col(b))
        .then(nw.col(a) + "_" + nw.col(b))
        .otherwise(nw.col(b) + "_" + nw.col(a))
    )
