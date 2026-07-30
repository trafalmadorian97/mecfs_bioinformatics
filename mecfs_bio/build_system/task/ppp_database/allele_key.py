"""
Order-agnostic allele key, shared across the PPP-database tasks for allele-aware
variant matching (so a swapped effect/non-effect orientation still matches).
"""

import polars as pl


def unordered_allele_key(a: str, b: str) -> pl.Expr:
    """Sort the two allele columns and join them, so {A, B} == {B, A}."""
    return (
        pl.when(pl.col(a) <= pl.col(b))
        .then(pl.col(a) + pl.lit("_") + pl.col(b))
        .otherwise(pl.col(b) + pl.lit("_") + pl.col(a))
    )
