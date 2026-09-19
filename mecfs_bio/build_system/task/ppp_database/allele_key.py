"""
Order-agnostic allele key, shared across the PPP-database tasks for allele-aware
variant matching (so a swapped effect/non-effect orientation still matches).

The key is valid only for SNVs. Mirrored indels such as T/TCA and TCA/T are distinct
variants, yet they sort to the same {A, B} key and would collide. Callers that key on it
must guarantee SNV-only input.

"""

import polars as pl


def unordered_allele_key(a: str, b: str) -> pl.Expr:
    """Sort the two allele columns and join them, so {A, B} == {B, A}.

    Valid only for SNVs: mirrored indels (T/TCA vs TCA/T) sort to the same key though they
    are distinct variants. Guard callers with assert_all_snv."""
    return (
        pl.when(pl.col(a) <= pl.col(b))
        .then(pl.col(a) + pl.lit("_") + pl.col(b))
        .otherwise(pl.col(b) + pl.lit("_") + pl.col(a))
    )


def assert_all_snv(df: pl.DataFrame, *allele_cols: str) -> None:
    """Fail fast if any allele in the named columns is not a single base."""
    assert allele_cols, "assert_all_snv requires at least one allele column"
    non_snv = df.filter(
        pl.any_horizontal(pl.col(c).str.len_bytes() != 1 for c in allele_cols)
    )
    assert non_snv.height == 0, (
        f"{non_snv.height} non-SNV variant(s) in columns {allele_cols}; the unordered "
        f"allele key is invalid for indels. First rows:\n{non_snv.select(allele_cols).head()}"
    )
